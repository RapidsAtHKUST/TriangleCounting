//
// Created by yche on 10/10/19.
//

#include "util/file_system/file_util.h"
#include "util/log.h"
#include "util/program_options/popl.h"
#include "util/search_util.h"

#include <chrono>
#include <cassert>

#define ATOMIC

using namespace std;
using namespace popl;
using namespace std::chrono;

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

int main(int argc, char *argv[]) {
    OptionParser op("Allowed options");
    auto string_option = op.add<Value<std::string>>("f", "file-path", "the graph bin file path");

    op.parse(argc, argv);
    using Edge = pair<int32_t, int32_t>;
    auto load_start = high_resolution_clock::now();

    if (string_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);

        auto file_name = string_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
        Edge *edge_lst =
                (Edge *) mmap(0, size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, file_fd, 0);

        int32_t max_node_id = -1;
        vector<int32_t> deg_lst;
        vector<uint32_t> off;
        vector<int32_t> adj_lst;
        vector<uint32_t> cur_write_off;
        auto load_end = high_resolution_clock::now();
        log_info("load edge list bin time: %.3lf s",
                 duration_cast<milliseconds>(load_end - load_start).count() / 1000.0);

        auto start = high_resolution_clock::now();
#if defined(LOCKS)
        omp_lock_t *locks;
#endif

#pragma omp parallel
        {
            // 1st: get the cardinality of degree array
#pragma omp for reduction(max: max_node_id)
            for (uint32_t i = 0; i < num_edges; i++) {
                max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
            }
#pragma omp single nowait
            {
                deg_lst = vector<int32_t>(static_cast<uint32_t>(max_node_id + 1));

            }
#pragma omp single nowait
            {
                off = vector<uint32_t>(static_cast<uint32_t>(max_node_id + 2));
                off[0] = 0;
            }
#pragma omp for
            for (auto i = 0; i < deg_lst.size(); i++) {
                deg_lst[i] = 0;
            }

            // 2nd: to count grouped neighbors, store in `deg_lst`
#pragma omp for
            for (uint32_t i = 0; i < num_edges; i++) {
                // atomic add for edge.first
                auto src = edge_lst[i].first;
                auto dst = edge_lst[i].second;
                int inc_deg_val, cur_deg_src;
                do {
                    cur_deg_src = deg_lst[src];
                    inc_deg_val = cur_deg_src + 1;
                } while (!__sync_bool_compare_and_swap(&(deg_lst[src]), cur_deg_src, inc_deg_val));
                do {
                    cur_deg_src = deg_lst[dst];
                    inc_deg_val = cur_deg_src + 1;
                } while (!__sync_bool_compare_and_swap(&(deg_lst[dst]), cur_deg_src, inc_deg_val));
            }

            // 3rd: compute prefix_sum and then scatter
#pragma omp single
            {
                for (int i = 0; i < deg_lst.size(); i++) {
                    off[i + 1] = off[i] + deg_lst[i];
                }
            }

#pragma omp single nowait
            {
                adj_lst = vector<int32_t>(off[off.size() - 1]);
            }

#pragma omp single nowait
            {
                cur_write_off = off;
            }

#pragma omp single nowait
            {
#if defined(LOCKS)
                locks = new omp_lock_t[deg_lst.size()];
#endif
            }

#if defined(LOCKS)
#pragma omp for
            for (int i = 0; i < deg_lst.size(); i++) {
                omp_init_lock(&locks[i]);
            }
#endif

            // 4th: barrier before we do the computation, and then construct destination vertices in CSR
#pragma omp single
            {
                auto middle = high_resolution_clock::now();
                log_info("before csr transform time: %.3lf s",
                         duration_cast<milliseconds>(middle - start).count() / 1000.0);
            }

#pragma omp for
            for (uint32_t i = 0; i < num_edges; i++) {
                auto src = edge_lst[i].first;
                auto dst = edge_lst[i].second;

                uint32_t new_offset, old_offset;
#if defined(ATOMIC)
                do {
                    old_offset = cur_write_off[src];
                    new_offset = old_offset + 1;
                } while (!__sync_bool_compare_and_swap(&(cur_write_off[src]), old_offset, new_offset));
                adj_lst[old_offset] = dst;

                do {
                    old_offset = cur_write_off[dst];
                    new_offset = old_offset + 1;
                } while (!__sync_bool_compare_and_swap(&(cur_write_off[dst]), old_offset, new_offset));
                adj_lst[old_offset] = src;
#elif defined(LOCKS)
                omp_set_lock(&locks[src]);
            old_offset = cur_write_off[src];
            new_offset = old_offset + 1;
            cur_write_off[src] = new_offset;
            omp_unset_lock(&locks[src]);

            adj_lst[old_offset] = dst;

            omp_set_lock(&locks[dst]);
            old_offset = cur_write_off[dst];
            new_offset = old_offset + 1;
            cur_write_off[dst] = new_offset;
            omp_unset_lock(&locks[dst]);

            adj_lst[old_offset] = src;
#endif
            }
#pragma omp single
            {
                auto middle2 = high_resolution_clock::now();
                log_info("before sort time: %.3lf s",
                         duration_cast<milliseconds>(middle2 - start).count() / 1000.0);
            }
            // 5th: sort each local ranges
#pragma omp for schedule(dynamic)
            for (auto u = 0; u < deg_lst.size(); u++) {
                assert(cur_write_off[u] == off[u + 1]);
                sort(begin(adj_lst) + off[u], begin(adj_lst) + off[u + 1]);
            }
        }
        auto end = high_resolution_clock::now();
        log_info("edge list to csr time: %.3lf s", duration_cast<milliseconds>(end - start).count() / 1000.0);

//        vector<uint32_t> off;
//        vector<int32_t> adj_lst;
        graph_t g;
        g.n = off.size()-1;
        g.m = adj_lst.size();
        g.adj = &adj_lst.front();
        g.num_edges = &off.front();
        size_t tc_cnt =0;
#pragma omp parallel for schedule(dynamic, 6000) reduction(+:tc_cnt)
        for (auto i = 0u; i < g.m; i++)
            ComputeSupport(&g, tc_cnt, i);
        log_info("There are %zu triangles in the input graph.", tc_cnt);
    }
}