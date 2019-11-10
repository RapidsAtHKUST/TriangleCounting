//
// Created by yche on 10/10/19.
//

#include <chrono>
#include <cassert>

#include <omp.h>

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

        // Remove Multi-Edges and Self-Loops.
        log_info("[%s] Mem: %d KB", __FUNCTION__, getValue());
        auto max_node_id = RemoveDuplicates(edge_lst, num_edges, edge_lst_buffer);
        log_info("[%s] Mem: %d KB", __FUNCTION__, getValue());

        auto num_vertices = static_cast<uint32_t >(max_node_id) + 1;
        log_info("Load Edge List Time: %.9lf s", global_timer.elapsed());

        // Convert Edge List to CSR.
        graph_t g{.n=num_vertices, .m = static_cast<long>(2L * num_edges),
                .adj=nullptr, .row_ptrs=nullptr};
        uint32_t *deg_lst;
        log_info("Undirected Graph G = (|V|, |E|): %lld, %lld", g.n, g.m / 2);

        g.adj = reinterpret_cast<int32_t *>(edge_lst_buffer);
#ifdef DODG
        ConvertEdgeListToDODGCSR(num_edges, edge_lst, num_vertices, deg_lst, g.row_ptrs, g.adj, max_omp_threads);
        assert(g.row_ptrs[num_vertices] == num_edges);
#else
        ConvertEdgeListToCSR(num_edges, edge_lst, num_vertices, deg_lst, g.row_ptrs, g.adj, max_omp_threads);
        assert(g.row_ptrs[num_vertices] == 2 * num_edges);
        log_debug("%d, %d", g.row_ptrs[num_vertices], 2 * num_edges);
#endif

        vector<int32_t> new_dict;
        vector<int32_t> old_dict;
        auto *tmp_mem_blocks = reinterpret_cast<int32_t *>(edge_lst);
#ifdef DODG
        ReorderDegDescendingDODG(g, new_dict, old_dict, tmp_mem_blocks, deg_lst);
#else
        ReorderDegDescending(g, new_dict, old_dict, tmp_mem_blocks);
#endif
        log_info("[%s] Mem: %d KB", __FUNCTION__, getValue());
        free(tmp_mem_blocks);
        log_info("[%s] Mem: %d KB", __FUNCTION__, getValue());

        // All-Edge Triangle Counting.
        size_t tc_cnt = 0;
        tc_cnt = CountTriBMPAndMergeWithPack(g, max_omp_threads);

        log_info("There are %zu triangles in the input graph.", tc_cnt);
        printf("There are %zu triangles in the input graph.\n", tc_cnt);
    }
}