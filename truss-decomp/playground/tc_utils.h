#pragma once

#include <vector>

#include <util/search_util.h>

using namespace std;

typedef struct {
    long n;
    long m;

    int32_t *adj;
    uint32_t *num_edges;
} graph_t;

inline int FindSrc(graph_t *g, int u, uint32_t edge_idx) {
    if (edge_idx >= g->num_edges[u + 1]) {
        // update last_u, preferring galloping instead of binary search because not large range here
        u = GallopingSearch(g->num_edges, static_cast<uint32_t>(u) + 1, g->n + 1, edge_idx);
        // 1) first > , 2) has neighbor
        if (g->num_edges[u] > edge_idx) {
            while (g->num_edges[u] - g->num_edges[u - 1] == 0) { u--; }
            u--;
        } else {
            // g->num_edges[u] == i
            while (g->num_edges[u + 1] - g->num_edges[u] == 0) {
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
            for (auto offset = g->num_edges[last_u]; offset < g->num_edges[last_u + 1]; offset++) {
                bits_vec[g->adj[offset]] = false;
            }
        }
        for (auto offset = g->num_edges[u]; offset < g->num_edges[u + 1]; offset++) {
            bits_vec[g->adj[offset]] = true;
        }
        last_u = u;
    }
    auto v = g->adj[i];
    auto du = g->num_edges[u + 1] - g->num_edges[u];
    auto dv = g->num_edges[v + 1] - g->num_edges[v];

    auto cnt = 0;
    if (du > dv || ((du == dv) && (u < v))) {
        cnt = ComputeCNHashBitVec(g, g->num_edges[v], g->num_edges[v + 1], bits_vec);
        tc_cnt += cnt;
    }
    return cnt;
}

// ru < rv, u points to v
inline bool less_than(int u, int v, int du, int dv) {
    return du > dv || ((du == dv) && (u < v));
}
