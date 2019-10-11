//
// Created by yche on 10/10/19.
//


#include <chrono>
#include <cassert>

#include <omp.h>

#include "util/file_system/file_util.h"
#include "util/log.h"
#include "util/program_options/popl.h"
#include "util/sort/parasort_cmp.h"
#include "util/timer.h"

#include "tc_utils.h"
#include "primitives.h"

using namespace std;
using namespace popl;
using namespace std::chrono;

template<typename T, typename I>
T RemoveDuplicates(pair<T, T> *&edge_lst, I &num_edges) {
    using Edge = pair<T, T>;
    Timer timer;
    T max_node_id = 0;
    // Partition.
#pragma omp parallel
    {
#pragma omp for reduction(max: max_node_id)
        for (auto i = 0u; i < num_edges; i++) {
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
        for (auto i = 0u; i < num_edges; i++) {
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
                (Edge *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_POPULATE, file_fd, 0);
        auto max_node_id = RemoveDuplicates(edge_lst, num_edges);

        auto deg_lst = vector<int32_t>(static_cast<uint32_t>(max_node_id + 1));
        auto off = vector<uint32_t>(static_cast<uint32_t>(max_node_id + 2));;
        vector<int32_t> adj_lst;
        vector<uint32_t> cur_write_off;
        auto load_end = high_resolution_clock::now();
        log_info("load edge list bin time: %.3lf s",
                 duration_cast<milliseconds>(load_end - load_start).count() / 1000.0);

        auto start = high_resolution_clock::now();

#pragma omp parallel
        {
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
            // 3rd: compute relative_off and then scatter
#pragma omp single
            {
                for (int i = 0; i < deg_lst.size(); i++) {
                    off[i + 1] = off[i] + deg_lst[i];
                }
                cur_write_off = off;
            }
#pragma omp single 
            {
                adj_lst = vector<int32_t>(off[off.size() - 1]);
            }

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

        graph_t g;
        g.n = off.size() - 1;
        g.m = adj_lst.size();
        g.adj = &adj_lst.front();
        g.num_edges = &off.front();
        size_t tc_cnt = 0;
#pragma omp parallel for schedule(dynamic, 6000) reduction(+:tc_cnt)
        for (auto i = 0u; i < g.m; i++)
            ComputeSupport(&g, tc_cnt, i);
        log_info("graph :%lld, %lld", g.n, g.m);
        log_info("There are %zu triangles in the input graph.", tc_cnt / 3);
    }
}