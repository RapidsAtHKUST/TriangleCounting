#pragma once

#include "primitives.h"
#include "graph.h"

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
            }
            max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
        }
#pragma omp single
        {
            num_buckets = max_node_id + 1;
        }
        // Partition.
        BucketSort(histogram, edge_lst, edge_lst_buffer, cur_write_off, bucket_ptrs, num_edges, num_buckets,
                   [&edge_lst](int i) {
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
    auto *relative_off = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
#pragma omp parallel num_threads(max_omp_threads)
    {
        SelectNotFOMP(histogram, edge_lst_buffer, edge_lst, relative_off, num_edges, [edge_lst](uint32_t it) {
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
        MemSetOMP(deg_lst, 0, num_vertices + 1);
        MemSetOMP(off, 0, num_vertices + 1);
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


inline void Reorder(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict, int32_t *&new_adj) {
    Timer timer;

    new_vid_dict = vector<int32_t>(g.n);
    vector<uint32_t> new_off(g.n + 1);
    new_off[0] = 0;

    auto max_omp_threads = omp_get_max_threads();
    auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);

#pragma omp parallel num_threads(max_omp_threads)
    {
        // 1st CSR: new_off, new_adj
#pragma omp for
        for (auto i = 0; i < g.n; i++) {
            new_vid_dict[old_vid_dict[i]] = i;
        }
        InclusivePrefixSumOMP(histogram, &new_off.front() + 1, g.n, [&g, &old_vid_dict](uint32_t new_id) {
            auto vertex = old_vid_dict[new_id];
            return g.row_ptrs[vertex + 1] - g.row_ptrs[vertex];
        });
#pragma omp single
        log_info("[%s]: Finish PrefixSum Time: %.9lf s", __FUNCTION__, timer.elapsed_and_reset());

        // 2nd Parallel Transform
#pragma omp for schedule(dynamic, 100)
        for (auto i = 0; i < g.n; i++) {
            auto origin_i = old_vid_dict[i];
            // transform
            auto cur_idx = new_off[i];
            for (auto my_old_off = g.row_ptrs[origin_i]; my_old_off < g.row_ptrs[origin_i + 1]; my_old_off++) {
                new_adj[cur_idx] = new_vid_dict[g.adj[my_old_off]];
                cur_idx++;
            }
            // sort the local ranges
            sort(new_adj + new_off[i], new_adj + new_off[i + 1]);
        }

        MemCpyOMP(g.row_ptrs, &new_off.front(), (g.n + 1));
    }
    swap(g.adj, new_adj);
    log_info("[%s]: Finish Reorder Time: %.3lf s", __FUNCTION__, timer.elapsed());
}

inline void ReorderDegDescending(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict,
                                 int32_t *&new_adj) {
    Timer timer;

#define USE_BUCKET_SORT
#ifdef USE_BUCKET_SORT
    auto max_omp_threads = omp_get_max_threads();
    auto max_deg = 0;
    auto *old_vid_dict_buffer = (int32_t *) malloc(sizeof(int32_t) * g.n);
    uint32_t *write_off = nullptr;
    uint32_t *bucket_ptrs = nullptr;
    auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);

#pragma omp parallel num_threads(max_omp_threads)
    {
#pragma omp for reduction(max: max_deg)
        for (auto i = 0; i < g.n; i++) {
            max_deg = max<int>(max_deg, g.row_ptrs[i + 1] - g.row_ptrs[i]);
        }
#pragma omp single nowait
        {
            old_vid_dict = vector<int32_t>(g.n);
        }
#pragma omp for
        for (auto i = 0u; i < g.n; i++) {
            old_vid_dict_buffer[i] = i;
        }
        auto ptr = &old_vid_dict[0];
        BucketSortSmallBuckets(histogram, old_vid_dict_buffer, ptr, write_off, bucket_ptrs,
                               g.n, max_deg + 1, [&g, old_vid_dict_buffer, max_deg](int i) {
                    auto u = old_vid_dict_buffer[i];
                    return max_deg - (g.row_ptrs[u + 1] - g.row_ptrs[u]);
                });
    }
    free(write_off);
    free(bucket_ptrs);
    free(old_vid_dict_buffer);
#else
    log_info("Use parallel sort (parasort)");
    old_vid_dict = vector<int32_t>(g.n);
#pragma omp parallel for
    for (auto i = 0; i < g.n; i++) {
        old_vid_dict[i] = i;
    }
    log_info("Allocation time:  %.9lf s", timer.elapsed());
    parasort(old_vid_dict.size(), &old_vid_dict.front(),
             [&g](int l, int r) -> bool {
                 return g.row_ptrs[l + 1] - g.row_ptrs[l] > g.row_ptrs[r + 1] - g.row_ptrs[r];
             },
             omp_get_max_threads());
#endif
    log_info("Deg-descending time:  %.9lf s", timer.elapsed());

    Reorder(g, new_vid_dict, old_vid_dict, new_adj);
}