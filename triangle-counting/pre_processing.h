#pragma once

#include "primitives.h"

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
#ifdef NO_REORDERING
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
#endif
    }
    free(cur_write_off);
    log_info("edge list to csr time: %.9lf s", convert_timer.elapsed());
}