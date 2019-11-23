//
// Created by yche on 11/10/19.
//

#include <cassert>

#include "ips4o/ips4o.hpp"

#include "util/pretty_print.h"
#include "util/popl.h"
#include "util/timer.h"
#include "util/util.h"
#include "util/log.h"
#include "util/primitives.h"

using namespace std;
using namespace popl;
using namespace std::chrono;

int main(int argc, char *argv[]) {
    OptionParser op("Allowed options");
    auto string_option = op.add<Value<std::string>>("f", "file-path", "the graph bin file path");
    op.parse(argc, argv);

    using Edge = pair<int32_t, int32_t>;
    Timer global_timer;

    if (string_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);

        auto file_name = string_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
        Edge *edge_lst = (Edge *) malloc(size);
#ifndef USE_LOG
        omp_set_num_threads(std::thread::hardware_concurrency());
#endif
        auto max_omp_threads = omp_get_max_threads();
#pragma omp parallel num_threads(max_omp_threads)
        {
            auto tid = omp_get_thread_num();
#pragma omp single
            log_info("Populate Mem Time: %.9lfs", global_timer.elapsed());

            size_t task_num = size;
            size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
            auto it_beg = avg * tid;
            auto it_end = min(avg * (tid + 1), task_num);
            auto *chars = reinterpret_cast<uint8_t *>(edge_lst);
            pread(file_fd, chars + it_beg, it_end - it_beg, it_beg);
        }
        log_info("Load File Time: %.9lfs", global_timer.elapsed());

        // Parallel Sort Globally.
        Timer sort_timer;
        ips4o::parallel::sort(edge_lst, edge_lst + num_edges, [](Edge l, Edge r) {
            if (l.first == r.first) {
                return l.second < r.second;
            }
            return l.first < r.first;
        });
        log_info("Sort Time: %.9lfs", sort_timer.elapsed());

        // Verification.
        int32_t max_vid = 0;
        int32_t num_buckets = 0;
        uint32_t *deg_histogram = nullptr;
        uint32_t max_d_dg = 0;
        uint32_t non_duplicates = 0;
#pragma omp parallel
        {
#pragma omp for
            for (size_t i = 0; i < num_edges - 1; i++) {
                auto l = edge_lst[i];
                auto r = edge_lst[i + 1];
                if (l.first == r.first) {
                    // = because of duplicates: self-loop-edge and multi-edge
                    assert(l.second <= r.second);
                } else {
                    assert(l.first < r.first);
                }
            }
#pragma omp single
            log_info("Verification Time: %.9lfs", sort_timer.elapsed());
#pragma omp for reduction(max: max_vid)
            for (size_t i = 0; i < num_edges; i++) {
                max_vid = max(max_vid, max(edge_lst[i].first, edge_lst[i].second));
            }
#pragma omp single
            {
                num_buckets = max_vid + 1;
                deg_histogram = (uint32_t *) malloc(sizeof(uint32_t) * num_buckets);
            }
            MemSetOMP(deg_histogram, 0, num_buckets);

            // Operator Fusion. (Do not select for the memory saving)
            auto duplicate_filter = [edge_lst](size_t it) {
                return edge_lst[it].first == edge_lst[it].second
                       || (it > 0 && edge_lst[it - 1] == edge_lst[it]);
            };
#pragma omp for reduction(+:non_duplicates)
            for (size_t i = 0; i < num_edges; i++) {
                if (!duplicate_filter(i)) {
                    non_duplicates++;
                    __sync_fetch_and_add(deg_histogram + edge_lst[i].first, 1);
                    __sync_fetch_and_add(deg_histogram + edge_lst[i].second, 1);
                }
            }
#pragma omp for reduction(max:max_d_dg)
            for (size_t i = 0; i < num_buckets; i++) {
                max_d_dg = max(max_d_dg, deg_histogram[i]);
            }
        }

        // Deg-Descending Dictionary.
        log_info("max-d: %d", max_d_dg);
        log_info("Use parallel sort (parasort)");
        auto old_vid_dict = vector<int32_t>(num_buckets);
#pragma omp parallel for
        for (auto i = 0; i < num_buckets; i++) {
            old_vid_dict[i] = i;
        }
        ips4o::parallel::sort(begin(old_vid_dict), end(old_vid_dict),
                              [deg_histogram](int l, int r) -> bool {
                                  return deg_histogram[l] > deg_histogram[r];
                              },
                              omp_get_max_threads());
        size_t num_large_deg_edges = 0;
        for (size_t i = 0; i < 32768; i++) {
            num_large_deg_edges += deg_histogram[old_vid_dict[i]];
            if (i % 1024 == 1024 - 1)
                log_info("%d, Large Deg Edges: %zu/%zu/%zu, Ratio: %.3lf", i,
                         num_large_deg_edges, non_duplicates, num_edges,
                         static_cast<double>(num_large_deg_edges) / non_duplicates);
        }
    }
}