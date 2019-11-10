#pragma once

#include "pre_processing.h"

bool RankLT(int du, int dv, int u, int v) {
    assert(u != v);
    return du < dv || ((du == dv) && (u < v));
}

template<typename T>
void ConvertEdgeListToDODGCSR(uint32_t num_edges, pair<T, T> *edge_lst,
                              uint32_t num_vertices, uint32_t *&deg_lst, uint32_t *&off, int32_t *&adj_lst,
                              int max_omp_threads) {
    Timer convert_timer;
    deg_lst = (uint32_t *) (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    auto *dodg_deg_lst = (uint32_t *) (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    off = (uint32_t *) (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    auto cur_write_off = (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
    vector<uint32_t> histogram;
#pragma omp parallel num_threads(max_omp_threads)
    {
        MemSetOMP(deg_lst, 0, num_vertices + 1);
        MemSetOMP(off, 0, num_vertices + 1);
#pragma omp single
        log_info("[%s]: InitTime: %.9lf s", __FUNCTION__, convert_timer.elapsed());

        // Histogram.
        EdgeListHistogram(num_vertices, num_edges, edge_lst, deg_lst, &convert_timer);
        MemSetOMP(dodg_deg_lst, 0, num_vertices + 1);
#pragma omp for
        for (uint32_t i = 0u; i < num_edges; i++) {
            auto src = edge_lst[i].first;
            auto dst = edge_lst[i].second;
            if (RankLT(deg_lst[src], deg_lst[dst], src, dst))
                __sync_fetch_and_add(&dodg_deg_lst[src], 1);
            else
                __sync_fetch_and_add(&dodg_deg_lst[dst], 1);
        }

#pragma omp single
        log_info("[%s]: Histogram Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());

        // PrefixSum.
        InclusivePrefixSumOMP(histogram, off + 1, num_vertices, [&dodg_deg_lst](uint32_t it) {
            return dodg_deg_lst[it];
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
            if (!RankLT(deg_lst[src], deg_lst[dst], src, dst)) {
                swap(src, dst);
            }
            auto old_offset = __sync_fetch_and_add(&(cur_write_off[src]), 1);
            adj_lst[old_offset] = dst;
        }
    }
    free(dodg_deg_lst);
    free(cur_write_off);
    log_info("[%s]: Total Conversion Time: %.9lf s", __FUNCTION__, convert_timer.elapsed());
}

inline void ReorderDegDescendingDODG(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict,
                                     int32_t *&new_adj, uint32_t *&deg_lst) {
    Timer timer;

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
            max_deg = max<int>(max_deg, deg_lst[i]);
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
                               g.n, max_deg + 1, [deg_lst, old_vid_dict_buffer, max_deg](int i) {
                    auto u = old_vid_dict_buffer[i];
                    return max_deg - (deg_lst[u]);
                });
    }
    free(write_off);
    free(bucket_ptrs);
    free(old_vid_dict_buffer);

    log_info("Deg-descending time:  %.9lf s", timer.elapsed());

    Reorder(g, new_vid_dict, old_vid_dict, new_adj);
}
