#pragma once

#include <vector>

#include "util/search_util.h"
#include "primitives.h"

using namespace std;

typedef struct {
    long n;
    long m;

    int32_t *adj;
    uint32_t *row_ptrs;
} graph_t;

// Assuming (offset_beg != offset_end)
template<typename T>
uint32_t GallopingSearchGeneral(T *array, uint32_t offset_beg, uint32_t offset_end, T val) {
    if (offset_beg == offset_end) { return offset_end; }
    if (array[offset_end - 1] < val) {
        return offset_end;
    }
    // galloping
    if (array[offset_beg] >= val) {
        return offset_beg;
    }
    if (array[offset_beg + 1] >= val) {
        return offset_beg + 1;
    }
    if (array[offset_beg + 2] >= val) {
        return offset_beg + 2;
    }

    auto jump_idx = 4u;
    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BranchFreeBinarySearch(array, (jump_idx >> 1u) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1u;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BranchFreeBinarySearch(array, (jump_idx >> 1u) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}

inline int FindSrc(graph_t *g, int u, uint32_t edge_idx) {
    if (edge_idx >= g->row_ptrs[u + 1]) {
        // update last_u, preferring galloping instead of binary search because not large range here
        u = GallopingSearchGeneral(g->row_ptrs, static_cast<uint32_t>(u) + 1, g->n + 1, edge_idx);
        // 1) first > , 2) has neighbor
        if (g->row_ptrs[u] > edge_idx) {
            while (g->row_ptrs[u] - g->row_ptrs[u - 1] == 0) { u--; }
            u--;
        } else {
            // g->row_ptrs[u] == i
            while (g->row_ptrs[u + 1] - g->row_ptrs[u] == 0) {
                u++;
            }
        }
    }
    return u;
}

template<typename T>
int ComputeCNHashBitVec(graph_t *g, uint32_t offset_beg, uint32_t offset_end, T &neighbor_bits) {
    auto cn_count = 0;
    for (auto offset = offset_beg; offset < offset_end; offset++) {
        if (neighbor_bits[g->adj[offset]]) {
            cn_count++;
        }
    }
    return cn_count;
}

inline int ComputeSupport(graph_t *g, size_t &tc_cnt, uint32_t i) {
    static thread_local auto u = 0;
    u = FindSrc(g, u, i);
    static thread_local auto last_u = -1;
    static thread_local auto bits_vec = vector<bool>(g->n, false);

    if (last_u != u) {
        // clear previous
        if (last_u != -1) {
            for (auto offset = g->row_ptrs[last_u]; offset < g->row_ptrs[last_u + 1]; offset++) {
                bits_vec[g->adj[offset]] = false;
            }
        }
        for (auto offset = g->row_ptrs[u]; offset < g->row_ptrs[u + 1]; offset++) {
            bits_vec[g->adj[offset]] = true;
        }
        last_u = u;
    }
    auto v = g->adj[i];
    auto du = g->row_ptrs[u + 1] - g->row_ptrs[u];
    auto dv = g->row_ptrs[v + 1] - g->row_ptrs[v];

    auto cnt = 0;
    if (du > dv || ((du == dv) && (u < v))) {
        cnt = ComputeCNHashBitVec(g, g->row_ptrs[v], g->row_ptrs[v + 1], bits_vec);
        tc_cnt += cnt;
    }
    return cnt;
}

void Reorder(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict, int32_t *&new_adj) {
    Timer timer;

    new_vid_dict = vector<int32_t>(g.n);
    vector<int32_t> new_deg(g.n);
    vector<uint32_t> new_off(g.n + 1);
    new_off[0] = 0;

    auto max_omp_threads = omp_get_max_threads();
    auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);

#pragma omp parallel num_threads(max_omp_threads)
    {
        auto tid = omp_get_thread_num();
        // 1st CSR: new_off, new_adj
#pragma omp for
        for (auto i = 0; i < g.n; i++) {
            new_vid_dict[old_vid_dict[i]] = i;
        }
#pragma omp for
        for (auto new_id = 0; new_id < g.n; new_id++) {
            auto vertex = old_vid_dict[new_id];
            new_deg[new_id] = g.row_ptrs[vertex + 1] - g.row_ptrs[vertex];
        }
        InclusivePrefixSumOMP(histogram, &new_off.front() + 1, g.n, [&new_deg](uint32_t it) {
            return new_deg[it];
        }, max_omp_threads);
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

        MemCpyOMP(g.row_ptrs, &new_off.front(), (g.n + 1), tid, max_omp_threads);
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
//                    assert(u < g.n);
                    return max_deg - (g.row_ptrs[u + 1] - g.row_ptrs[u]);
                }, max_omp_threads);
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

