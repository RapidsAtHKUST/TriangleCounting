#pragma once

#include "util/primitives.h"
#include "ips4o/ips4o.hpp"
#include "timer.h"
#include "log.h"

template<typename T, typename OFF>
T RemoveDuplicates(pair<T, T> *&edge_lst, OFF &num_edges, pair<T, T> *&edge_lst_buffer) {
    using Edge = pair<T, T>;
    Timer timer;
    T max_node_id = 0;
    T num_buckets;
    auto max_omp_threads = omp_get_max_threads();

    OFF *bucket_ptrs;
    OFF *cur_write_off;
    vector<OFF> histogram;

#pragma omp parallel num_threads(max_omp_threads)
    {
#pragma omp for reduction(max: max_node_id)
        for (OFF i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
            }
            max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
        }
#pragma omp single
        {
            num_buckets = max_node_id + 1;
        }
        // Partition.
        BucketSort(histogram, edge_lst, edge_lst_buffer, cur_write_off, bucket_ptrs, num_edges, num_buckets,
                   [&edge_lst](size_t i) {
                       return edge_lst[i].first;
                   }, &timer);

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
    auto *relative_off = (OFF *) malloc(sizeof(OFF) * num_edges);
#pragma omp parallel num_threads(max_omp_threads)
    {
        SelectNotFOMP(histogram, edge_lst_buffer, edge_lst, relative_off, num_edges, [edge_lst](size_t it) {
            return edge_lst[it].first == edge_lst[it].second || (it > 0 && edge_lst[it - 1] == edge_lst[it]);
        });
    }
    swap(edge_lst, edge_lst_buffer);
    num_edges = num_edges - relative_off[num_edges - 1];
    free(relative_off);
    log_info("New # of edges: %zu, Elapsed: %.9lfs", num_edges, timer.elapsed());
    log_debug("max_node_id: %d", max_node_id);
    return max_node_id;
}

template<typename T, typename D, typename I, typename OFF, typename F>
void EdgeListHistogram(I num_vertices, OFF num_edges, pair<T, T> *edge_lst, D *deg_lst, F f) {
    auto local_buf = (uint8_t *) calloc(num_vertices, sizeof(uint8_t));
#pragma omp for
    for (size_t i = 0u; i < num_edges; i++) {
        if (f(i)) {
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
    }
    for (size_t i = 0u; i < num_vertices; i++) {
        // atomic add for edge.first
        if (local_buf[i] > 0)
            __sync_fetch_and_add(&(deg_lst[i]), local_buf[i]);
    }
    free(local_buf);
#pragma omp barrier
}

template<typename T, typename D, typename I, typename OFF>
void EdgeListHistogram(I num_vertices, OFF num_edges, pair<T, T> *edge_lst, D *deg_lst) {
    EdgeListHistogram(num_vertices, num_edges, edge_lst, deg_lst, [](size_t it) {
        return true;
    });
}


template<typename T, typename OFF>
void ConvertEdgeListToCSR(OFF num_edges, pair<T, T> *edge_lst,
                          uint32_t num_vertices, uint32_t *&deg_lst, OFF *&off,
                          int32_t *&adj_lst, int max_omp_threads) {
    Timer convert_timer;
    deg_lst = (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    off = (OFF *) malloc(sizeof(OFF) * (num_vertices + 1));
    auto cur_write_off = (OFF *) malloc(sizeof(OFF) * (num_vertices + 1));
    vector<OFF> histogram;
#pragma omp parallel num_threads(max_omp_threads)
    {
        MemSetOMP(deg_lst, 0, num_vertices + 1);
        MemSetOMP(off, 0, num_vertices + 1);
#pragma omp single
        log_info("[%s]: InitTime: %.9lf s", __FUNCTION__, convert_timer.elapsed());
        EdgeListHistogram(num_vertices, num_edges, edge_lst, deg_lst);

#pragma omp single
        log_info("[%s]: Histogram Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());

        // PrefixSum.
        InclusivePrefixSumOMP(histogram, off + 1, num_vertices, [&deg_lst](uint32_t it) {
            return deg_lst[it];
        });
        MemCpyOMP(cur_write_off, off, num_vertices + 1);

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
        for (size_t i = 0; i < num_edges; i++) {
            auto src = edge_lst[i].first;
            auto dst = edge_lst[i].second;
            auto old_offset = __sync_fetch_and_add(&(cur_write_off[src]), 1);
            adj_lst[old_offset] = dst;
            old_offset = __sync_fetch_and_add(&(cur_write_off[dst]), 1);
            adj_lst[old_offset] = src;
        }
    }
    free(cur_write_off);
    log_info("[%s]: Total Conversion Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());
}
