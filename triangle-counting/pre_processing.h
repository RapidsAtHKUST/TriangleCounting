#pragma once

#include "util/primitives.h"
#include "util/graph.h"

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
#pragma omp barrier
}

inline void Reorder(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict, int32_t *&new_adj) {
    Timer timer;

    new_vid_dict = vector<int32_t>(g.n);
    vector<row_ptr_t> new_off(g.n + 1);
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