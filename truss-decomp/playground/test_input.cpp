//
// Created by yche on 10/10/19.
//


#include <chrono>
#include <cassert>

#include <omp.h>

#include "util/file_system/file_util.h"
#include "util/log.h"
#include "util/program_options/popl.h"
#include "util/sort/parasort_cmp.h"
#include "util/timer.h"
#include "util/stat.h"
#include "util/util.h"

#include "tc_utils.h"
#include "primitives.h"

using namespace std;
using namespace popl;
using namespace std::chrono;

template<typename T, typename I>
T RemoveDuplicates(pair<T, T> *&edge_lst, I &num_edges) {
    using Edge = pair<T, T>;
    Timer timer;
    T max_node_id = 0;
    // Partition.
#pragma omp parallel
    {
#pragma omp for reduction(max: max_node_id)
        for (I i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
                max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
            }
        }
    }
    auto max_omp_threads = omp_get_max_threads();
    parasort(num_edges, edge_lst, [](const Edge &left, const Edge &right) {
        if (left.first == right.first) {
            return left.second < right.second;
        }
        return left.first < right.first;
    }, max_omp_threads);
    log_info("Finish Sort, %.9lfs", timer.elapsed());

    // Selection.
    auto *relative_off = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
    Edge *edge_lst2 = (Edge *) malloc(sizeof(Edge) * num_edges);
    auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);
#pragma omp parallel
    {
        FlagPrefixSumOMP(histogram, relative_off, num_edges, [edge_lst](uint32_t it) {
            return edge_lst[it].first == edge_lst[it].second || (it > 0 && edge_lst[it - 1] == edge_lst[it]);
        }, max_omp_threads);
#pragma omp for
        for (I i = 0u; i < num_edges; i++) {
            if (!(edge_lst[i].first == edge_lst[i].second || (i > 0 && edge_lst[i - 1] == edge_lst[i]))) {
                auto off = i - relative_off[i];
                edge_lst2[off] = edge_lst[i];
            }
        }
    }
    edge_lst = edge_lst2;
    num_edges = num_edges - relative_off[num_edges - 1];

    log_info("New # of edges: %zu, Elapsed: %.9lfs", num_edges, timer.elapsed());
    return max_node_id;
}

template<typename T>
void ConvertEdgeListToCSR(uint32_t num_edges, pair<T, T> *edge_lst,
                          uint32_t num_vertices, uint32_t *&deg_lst, uint32_t *&off,
                          int32_t *&adj_lst, int max_omp_threads) {
    Timer convert_timer;
    deg_lst = (uint32_t *) (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    off = (uint32_t *) (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    auto cur_write_off = (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);

#pragma omp parallel num_threads(max_omp_threads)
    {
        auto tid = omp_get_thread_num();
        {
            size_t task_num = num_vertices + 1;
            size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
            auto it_beg = avg * tid;
            auto it_end = min(avg * (tid + 1), task_num);
            memset(deg_lst + it_beg, 0, sizeof(uint32_t) * (it_end - it_beg));
            memset(off + it_beg, 0, sizeof(uint32_t) * (it_end - it_beg));
#pragma omp barrier
        }
        // Histogram.
#pragma omp for
        for (uint32_t i = 0u; i < num_edges; i++) {
            // atomic add for edge.first
            auto src = edge_lst[i].first;
            auto dst = edge_lst[i].second;
            __sync_fetch_and_add(&(deg_lst[src]), 1);
            __sync_fetch_and_add(&(deg_lst[dst]), 1);
        }
        // PrefixSum.
        {
            InclusivePrefixSumOMP(histogram, off + 1, num_vertices, [&deg_lst](uint32_t it) {
                return deg_lst[it];
            }, max_omp_threads);
            size_t task_num = num_vertices + 1;
            size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
            auto it_beg = avg * tid;
            auto it_end = min(avg * (tid + 1), task_num);
            memcpy(cur_write_off + it_beg, off + it_beg, sizeof(uint32_t) * (it_end - it_beg));
#pragma omp barrier
        }
#pragma omp single
        {
            adj_lst = (int32_t *) malloc(sizeof(int32_t) * off[num_vertices]);
            log_info("before csr transform time: %.9lf s", convert_timer.elapsed());
        }
        // Scatter.
#pragma omp for
        for (uint32_t i = 0; i < num_edges; i++) {
            auto src = edge_lst[i].first;
            auto dst = edge_lst[i].second;
            auto old_offset = __sync_fetch_and_add(&(cur_write_off[src]), 1);
            adj_lst[old_offset] = dst;
            old_offset = __sync_fetch_and_add(&(cur_write_off[dst]), 1);
            adj_lst[old_offset] = src;
        }
#pragma omp single
        {
            log_info("before sort time: %.9lf s", convert_timer.elapsed());
        }
        // 5th: sort each local ranges
#pragma omp for schedule(dynamic)
        for (auto u = 0u; u < num_vertices; u++) {
            assert(cur_write_off[u] == off[u + 1]);
            sort(adj_lst + off[u], adj_lst + off[u + 1]);
        }
    }
    free(cur_write_off);
    log_info("edge list to csr time: %.9lf s", convert_timer.elapsed());
}

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
        Edge *edge_lst = (Edge *) mmap(nullptr, size, PROT_READ | PROT_WRITE,
                                       MAP_PRIVATE | MAP_POPULATE, file_fd, 0);

        // Remove Multi-Edges and Self-Loops.
//        Edge *prev_edge_lst = edge_lst;
//        auto prev_num_edges = num_edges;
        auto max_node_id = RemoveDuplicates(edge_lst, num_edges);
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
        assert(g.row_ptrs[num_vertices] == 2 * num_edges);

        vector<int32_t> new_dict;
        vector<int32_t> old_dict;
        ReorderDegDescending(g, new_dict, old_dict);

        // All-Edge Triangle Counting.
        size_t tc_cnt = 0;
        int max_d = 0;
#ifdef BASELINE
#pragma omp parallel for schedule(dynamic, 6000) reduction(+:tc_cnt)
        for (auto i = 0u; i < g.m; i++)
            ComputeSupport(&g, tc_cnt, i);
        tc_cnt /=3;
#else
        Timer tc_timer;

        auto *row_ptrs_end = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
#pragma omp parallel num_threads(max_omp_threads)
        {
            // Bit vectors.
            auto bits_vec = vector<bool>(g.n, false);
#pragma omp for reduction(max: max_d)
            for (auto u = 0u; u < g.n; u++) {
                row_ptrs_end[u + 1] = static_cast<uint32_t>(
                        lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1], u) - g.adj);
                max_d = max<int>(max_d, row_ptrs_end[u + 1] - g.row_ptrs[u]);
            }
#pragma omp single
            log_info("finish init row_ptrs_end, max d: %d", max_d);

#pragma omp for schedule(dynamic, 100) reduction(+:tc_cnt)
            for (auto u = 0u; u < g.n; u++) {
                // Set.
                for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                    auto v = g.adj[edge_idx];
                    bits_vec[v] = true;
                }
                for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                    auto v = g.adj[edge_idx];
                    tc_cnt += ComputeCNHashBitVec(&g, g.row_ptrs[v], row_ptrs_end[v + 1], bits_vec);
                }
                // Clear.
                for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                    auto v = g.adj[edge_idx];
                    bits_vec[v] = false;
                }
            }
        }
        free(row_ptrs_end);
        log_info("Forward cost: %.3lf s, Mem Usage: %s KB",
                 tc_timer.elapsed(), FormatWithCommas(getValue()).c_str());
        log_info("Triangle Cnt: %s", FormatWithCommas(tc_cnt).c_str());
#endif

        log_info("There are %zu triangles in the input graph.", tc_cnt);
    }
}