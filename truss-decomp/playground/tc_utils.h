#pragma once

#include <vector>

#include <util/search_util.h>

using namespace std;

typedef struct {
    long n;
    long m;

    int32_t *adj;
    uint32_t *row_ptrs;
} graph_t;

inline int FindSrc(graph_t *g, int u, uint32_t edge_idx) {
    if (edge_idx >= g->row_ptrs[u + 1]) {
        // update last_u, preferring galloping instead of binary search because not large range here
        u = GallopingSearch(g->row_ptrs, static_cast<uint32_t>(u) + 1, g->n + 1, edge_idx);
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

// ru < rv, u points to v
inline bool less_than(int u, int v, int du, int dv) {
    return du > dv || ((du == dv) && (u < v));
}

void Reorder(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict) {
    Timer timer;

    new_vid_dict = vector<int>(g.n);
    for (auto i = 0; i < g.n; i++) {
        new_vid_dict[old_vid_dict[i]] = i;
    }
    // new-deg
    vector<int> new_deg(g.n);
    for (auto new_id = 0; new_id < g.n; new_id++) {
        auto vertex = old_vid_dict[new_id];
        new_deg[new_id] = g.row_ptrs[vertex + 1] - g.row_ptrs[vertex];
        assert(new_deg[new_id] >= 0);
    }

    // verify permutation
    for (auto i = 0; i < std::min<int32_t>(5, static_cast<int32_t>(new_vid_dict.size())); i++) {
        log_info("old->new %d -> %d", i, new_vid_dict[i]);
    }
    vector<int> verify_map(new_vid_dict.size(), 0);
    int cnt = 0;

#pragma omp parallel
    {
#pragma omp for reduction(+:cnt)
        for (auto i = 0u; i < new_vid_dict.size(); i++) {
            if (verify_map[new_vid_dict[i]] == 0) {
                cnt++;
                verify_map[new_vid_dict[i]] = 1;
            } else {
                assert(false);
            }
        }
#pragma omp single
        log_info("%d, %d", cnt, new_vid_dict.size());
        assert(static_cast<size_t>(cnt) == new_vid_dict.size());
    }
    // 1st CSR: new_off, new_neighbors
    vector<uint32_t> new_off(g.n + 1);
    new_off[0] = 0;
    assert(static_cast<long>(new_off.size()) == g.n + 1);
    for (auto i = 0u; i < g.n; i++) { new_off[i + 1] = new_off[i] + new_deg[i]; }
    log_info("%zu", new_off[g.n]);
    assert(new_off[g.n] == g.m);

    vector<int> new_neighbors(g.m);

    log_info("init ordering structures time: %.9lf s", timer.elapsed_and_reset());

    // 2nd Parallel Transform
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 100)
        for (auto i = 0; i < g.n; i++) {
            auto origin_i = old_vid_dict[i];
            // transform
            auto cur_idx = new_off[i];
            for (auto my_old_off = g.row_ptrs[origin_i]; my_old_off < g.row_ptrs[origin_i + 1]; my_old_off++) {
                if (cur_idx > g.m) {
                    log_info("%d, i: %d", cur_idx, i);
                }
                assert(cur_idx <= g.m);
                assert(my_old_off <= g.m);
                assert(g.adj[my_old_off] < g.n);
                new_neighbors[cur_idx] = new_vid_dict[g.adj[my_old_off]];
                cur_idx++;
            }
            // sort the local ranges
            sort(begin(new_neighbors) + new_off[i], begin(new_neighbors) + new_off[i + 1]);
        }
    }
    log_info("parallel transform and sort: %.3lf s", timer.elapsed());

    memcpy(g.adj, &new_neighbors.front(), g.m * sizeof(int32_t));
    memcpy(g.row_ptrs, &new_off.front(), (g.n + 1) * sizeof(uint32_t));
    log_info("parallel transform after copying: %.3lf s", timer.elapsed());
}

inline void ReorderDegDescending(graph_t &g, vector<int32_t> &new_vid_dict, vector<int32_t> &old_vid_dict) {
    Timer timer;

    old_vid_dict = vector<int>(g.n);
    for (auto i = 0u; i < old_vid_dict.size(); i++) { old_vid_dict[i] = i; }

    log_info("Use parallel sort (parasort)");
    parasort(old_vid_dict.size(), &old_vid_dict.front(),
             [&g](int l, int r) -> bool {
                 return g.row_ptrs[l + 1] - g.row_ptrs[l] > g.row_ptrs[r + 1] - g.row_ptrs[r];
             },
             omp_get_max_threads());
    log_info("Deg-descending time:  %.9lf s", timer.elapsed());

    Reorder(g, new_vid_dict, old_vid_dict);
}

