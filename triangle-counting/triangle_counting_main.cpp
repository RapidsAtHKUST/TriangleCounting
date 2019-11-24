//
// Created by yche on 10/10/19.
//

#include <chrono>
#include <cassert>

#include <omp.h>
#include <malloc.h>

#include "ips4o/ips4o.hpp"

#include "util/log.h"
#include "util/popl.h"
#include "util/parasort_cmp.h"
#include "util/timer.h"
#include "util/util.h"
#include "util/pretty_print.h"

#include "pre_processing.h"
#include "pre_processing_dodg.h"
#include "triangle_counting.h"

using namespace std;
using namespace popl;
using namespace std::chrono;

//#define MMAP
#define PAGE_SIZE (4096)
#define IO_REQ_SIZE (PAGE_SIZE * 32)
#define IO_QUEUE_DEPTH (16)

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
#ifndef MMAP
        auto file_fd = open(file_name.c_str(), O_RDONLY | O_DIRECT, S_IRUSR | S_IWUSR);
        Edge *edge_lst = (Edge *) memalign(PAGE_SIZE, size + IO_REQ_SIZE);
#else
        auto file_fd = open(file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
        Edge *edge_lst = (Edge *) mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_PRIVATE, file_fd, 0);
        madvise(edge_lst, size, MADV_SEQUENTIAL);
#endif
#ifndef USE_LOG
        omp_set_num_threads(std::thread::hardware_concurrency());
#endif
        auto max_omp_threads = omp_get_max_threads();
#ifndef MMAP
        Timer io_timer;
        size_t read_size = 0;
#pragma omp parallel num_threads(IO_QUEUE_DEPTH)
        {
#pragma omp for schedule(dynamic, 1) reduction(+:read_size)
            for (size_t i = 0; i < size; i += IO_REQ_SIZE) {
                auto it_beg = i;
                auto *chars = reinterpret_cast<uint8_t *>(edge_lst);
                auto ret = pread(file_fd, chars + it_beg, IO_REQ_SIZE, it_beg);
                if (ret != IO_REQ_SIZE) {
                    log_error("Err, %zu, %zu, %zu, %d", i, it_beg, IO_REQ_SIZE, ret);
                } else {
                    read_size += ret;
                }
            }
#pragma omp single
            log_info("%zu, %zu", read_size, size);
        }
        log_info("IO Time: %.6lfs, DIO-QPS: %.6lf GB/s", io_timer.elapsed(), size / io_timer.elapsed() / pow(1024, 3));
#else
        mlock(edge_lst, size);
#endif
        log_info("Load File Time: %.9lfs", global_timer.elapsed());
        // 1st: Remove Multi-Edges and Self-Loops.
        Timer sort_timer;
        int32_t max_node_id = 0;
#pragma omp parallel for reduction(max: max_node_id) schedule(dynamic, 32*1024)
        for (size_t i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
            }
            max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
        }
        log_info("Populate File Time: %.9lfs", global_timer.elapsed());

        ips4o::parallel::sort(edge_lst, edge_lst + num_edges, [](Edge l, Edge r) {
            if (l.first == r.first) {
                return l.second < r.second;
            }
            return l.first < r.first;
        });

#ifdef VERIFY_PROCESSED_EDGE_LST
        size_t gt_num_edges = 0;
#pragma omp parallel
        {
#pragma omp for reduction(+:gt_num_edges)
            for (size_t it = 0; it < num_edges; it++) {
                if (edge_lst[it].first == edge_lst[it].second || (it > 0 && edge_lst[it - 1] == edge_lst[it]))
                    continue;
                gt_num_edges++;
            }
#pragma omp single
            {
                log_info("Remove Duplicates: %zu", gt_num_edges);
            }
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
        }
#endif
        log_info("Sort Time: %.9lfs", sort_timer.elapsed());
        auto num_vertices = static_cast<uint32_t >(max_node_id) + 1;
        log_info("Pre-Process Edge List Time: %.9lf s", global_timer.elapsed());

        // 2nd: Convert Edge List to CSR.
        graph_t g{.n=num_vertices, .m = 0, .adj=nullptr, .row_ptrs=nullptr};
        uint32_t *deg_lst;
        g.adj = (int32_t *) malloc(size / 2);

#ifdef MMAP
        munlock(edge_lst, size);
#endif
        ConvertEdgeListToDODGCSR(num_edges, edge_lst, num_vertices, deg_lst, g.row_ptrs, g.adj,
                                 max_omp_threads, [=](size_t it) {
                    return !(edge_lst[it].first == edge_lst[it].second
                             || (it > 0 && edge_lst[it - 1] == edge_lst[it]));
                });
        assert(g.row_ptrs[num_vertices] <= num_edges);
        g.m = g.row_ptrs[num_vertices];
        log_info("Undirected Graph G = (|V|, |E|): %lld, %lld", g.n, g.m);
        log_info("Mem Usage: %s KB", FormatWithCommas(getValue()).c_str());

        // 3rd: Reordering.
        vector<int32_t> new_dict;
        vector<int32_t> old_dict;
#ifndef MMAP
        free(edge_lst);
#else
        munmap(edge_lst, size);
#endif
        auto *tmp_mem_blocks = (int32_t *) malloc(size / 2);
        auto *org = g.adj;
        ReorderDegDescendingDODG(g, new_dict, old_dict, tmp_mem_blocks, deg_lst);
        free(org);

        // 4th: Triangle Counting.
        log_info("Mem Usage: %s KB", FormatWithCommas(getValue()).c_str());
        size_t tc_cnt = 0;
        tc_cnt = CountTriBMPAndMergeWithPackDODG(g, max_omp_threads);
        log_info("Mem Usage: %s KB", FormatWithCommas(getValue()).c_str());
        log_info("There are %zu triangles in the input graph.", tc_cnt);
        printf("There are %zu triangles in the input graph.\n", tc_cnt);
    }
}