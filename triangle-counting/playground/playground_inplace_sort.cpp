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
        Edge *edge_lst_buffer = (Edge *) malloc(size);
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

        Timer sort_timer;
        ips4o::parallel::sort(edge_lst, edge_lst + num_edges, [](Edge l, Edge r) {
            if (l.first == r.first) {
                return l.second < r.second;
            }
            return l.first < r.first;
        });
        log_info("Sort Time: %.9lfs", sort_timer.elapsed());
        // Verification
#pragma omp parallel
        {
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
        }
    }
}