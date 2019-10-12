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
#include "triangle_counting.h"

using namespace std;
using namespace popl;
using namespace std::chrono;


int main(int argc, char *argv[]) {
//#ifdef USE_LOG
//    setlocale(LC_NUMERIC, "");
//#endif
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
        Edge *edge_lst = (Edge *) mmap(nullptr, size, PROT_READ | PROT_WRITE,
                                       MAP_PRIVATE | MAP_POPULATE, file_fd, 0);

        // Remove Multi-Edges and Self-Loops.
        Edge *prev_edge_lst = edge_lst;
        auto max_node_id = RemoveDuplicates(edge_lst, num_edges);
        munmap(prev_edge_lst, size);

        auto num_vertices = static_cast<uint32_t >(max_node_id) + 1;
        log_info("load edge list bin time: %.9lf s", global_timer.elapsed());

        // Convert Edge List to CSR.
        graph_t g;
        g.n = num_vertices;
        g.m = 2L * num_edges;
        uint32_t *deg_lst;
        log_info("graph :%lld, %lld", g.n, g.m);

        auto max_omp_threads = omp_get_max_threads();
        ConvertEdgeListToCSR(num_edges, edge_lst, num_vertices, deg_lst, g.row_ptrs, g.adj, max_omp_threads);
        free(edge_lst);
        assert(g.row_ptrs[num_vertices] == 2 * num_edges);

        vector<int32_t> new_dict;
        vector<int32_t> old_dict;
        ReorderDegDescending(g, new_dict, old_dict);

        // All-Edge Triangle Counting.
        size_t tc_cnt = 0;
#ifdef BASELINE
#pragma omp parallel for schedule(dynamic, 6000) reduction(+:tc_cnt)
        for (auto i = 0u; i < g.m; i++)
            ComputeSupport(&g, tc_cnt, i);
        tc_cnt /=3;
#else
//        tc_cnt = CountTriBMP(g, max_omp_threads);
//        tc_cnt = CountTriBMPWithPack(g, max_omp_threads);
        tc_cnt = CountTriBMPAndMergeWithPack(g, max_omp_threads);
#endif

        log_info("There are %zu triangles in the input graph.", tc_cnt);
        printf("There are %zu triangles in the input graph.\n", tc_cnt);
    }
}