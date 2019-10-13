#pragma once

#include "primitives.h"

template<typename T, typename I>
T RemoveDuplicates(pair<T, T> *&edge_lst, I &num_edges, pair<T, T> *&edge_lst_buffer) {
    using Edge = pair<T, T>;
    Timer timer;
    T max_node_id = 0;
    T num_buckets;
    auto max_omp_threads = omp_get_max_threads();

    uint32_t *bucket_ptrs;
    uint32_t *cur_write_off;
    vector<uint32_t> histogram;

#pragma omp parallel num_threads(max_omp_threads)
    {
#pragma omp for reduction(max: max_node_id)
        for (I i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
                max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
            }
        }
#pragma omp single
        {
            num_buckets = max_node_id;
        }
        // Partition.
        BucketSort(histogram, edge_lst, edge_lst_buffer, cur_write_off, bucket_ptrs, num_edges, num_buckets,
                   [&edge_lst](int i) {
                       return edge_lst[i].first;
                   }, max_omp_threads, &timer);

        // Sort.
#pragma omp for schedule(dynamic, 600)
        for (auto i = 0; i < num_buckets; i++) {
            sort(edge_lst_buffer + bucket_ptrs[i], edge_lst_buffer + bucket_ptrs[i + 1],
                 [](const Edge &left, const Edge &right) {
                     return left.second < right.second;
                 });
        }
    }
    swap(edge_lst, edge_lst_buffer);
    free(cur_write_off);
    free(bucket_ptrs);
    log_info("Finish Sort, %.9lfs", timer.elapsed());

    // Selection.
    auto *relative_off = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
#pragma omp parallel num_threads(max_omp_threads)
    {
        SelectNotFOMP(histogram, edge_lst_buffer, edge_lst, relative_off, num_edges, [edge_lst](uint32_t it) {
            return edge_lst[it].first == edge_lst[it].second || (it > 0 && edge_lst[it - 1] == edge_lst[it]);
        }, max_omp_threads);
    }
    swap(edge_lst, edge_lst_buffer);
    num_edges = num_edges - relative_off[num_edges - 1];
    free(relative_off);
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
    vector<uint32_t> histogram;

#pragma omp parallel num_threads(max_omp_threads)
    {
        auto tid = omp_get_thread_num();
        MemSetOMP(deg_lst, 0, num_vertices + 1, tid, max_omp_threads);
        MemSetOMP(off, 0, num_vertices + 1, tid, max_omp_threads);
        auto local_buf = (uint8_t *) calloc(num_vertices, sizeof(uint8_t));
#pragma omp single
        log_info("[%s]: InitTime: %.9lf s", __FUNCTION__, convert_timer.elapsed());

        // Histogram.
#pragma omp for
        for (uint32_t i = 0u; i < num_edges; i++) {
            auto src = edge_lst[i].first;
            auto dst = edge_lst[i].second;
            local_buf[src]++;
            if (local_buf[src] == 0xff) {
                __sync_fetch_and_add(&deg_lst[src], 0xff);
                local_buf[src] = 0;
            }
            local_buf[dst]++;
            if (local_buf[dst] == 0xff) {
                __sync_fetch_and_add(&deg_lst[dst], 0xff);
                local_buf[dst] = 0;
            }
        }
#pragma omp single
        log_info("[%s]: Histogram Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());
        for (size_t i = 0; i < num_vertices; i++) {
            // atomic add for edge.first
            if (local_buf[i] > 0)
                __sync_fetch_and_add(&(deg_lst[i]), local_buf[i]);
        }
#pragma omp barrier

#pragma omp single
        log_info("[%s]: Histogram Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());

        // PrefixSum.
        InclusivePrefixSumOMP(histogram, off + 1, num_vertices, [&deg_lst](uint32_t it) {
            return deg_lst[it];
        }, max_omp_threads);
        MemCpyOMP(cur_write_off, off, num_vertices + 1, tid, max_omp_threads);

        // Scatter.
#pragma omp single
        {
            if (adj_lst == nullptr) {
                log_info("Allocate Inside (adj_lst)...");
                adj_lst = (int32_t *) malloc(sizeof(int32_t) * off[num_vertices]);
            }
            log_info("[%s]: PrefixSum Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());
        }

#pragma omp for
        for (uint32_t i = 0; i < num_edges; i++) {
            auto src = edge_lst[i].first;
            auto dst = edge_lst[i].second;
            auto old_offset = __sync_fetch_and_add(&(cur_write_off[src]), 1);
            adj_lst[old_offset] = dst;
            old_offset = __sync_fetch_and_add(&(cur_write_off[dst]), 1);
            adj_lst[old_offset] = src;
        }
        free(local_buf);
    }
    free(cur_write_off);
    log_info("[%s]: Total Conversion Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());
}