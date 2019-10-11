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
    Timer global_timer;
    if (string_option->is_set()) {
        size_t size = file_size(string_option->value(0).c_str());
        size_t num_edges = size / sizeof(uint32_t) / 2;
        log_info("File size: %zu", size);
        log_info("#of Edges: %zu", num_edges);

        auto file_name = string_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
        Edge *edge_lst =
                (Edge *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_POPULATE, file_fd, 0);

        // Remove Multi-Edges and Self-Loops.
        Edge *prev_edge_lst = edge_lst;
        auto prev_num_edges = num_edges;
        auto max_node_id = RemoveDuplicates(edge_lst, num_edges);

        auto max_omp_threads = omp_get_max_threads();
        auto histogram = vector<uint32_t>((max_omp_threads + 1) * CACHE_LINE_ENTRY, 0);
        auto num_vertices = max_node_id + 1;
        auto deg_lst = vector<int32_t>(static_cast<uint32_t>(num_vertices));
        auto off = vector<uint32_t>(static_cast<uint32_t>(num_vertices + 1));
        auto cur_write_off = (uint32_t *) malloc(sizeof(uint32_t) * (num_vertices + 1));
        int32_t *adj_lst;
        log_info("load edge list bin time: %.9lf s", global_timer.elapsed());

        Timer convert_timer;
#pragma omp parallel num_threads(max_omp_threads)
        {
            auto tid = omp_get_thread_num();
            // Histogram.
#pragma omp for
            for (uint32_t i = 0; i < num_edges; i++) {
                // atomic add for edge.first
                auto src = edge_lst[i].first;
                auto dst = edge_lst[i].second;
                __sync_fetch_and_add(&(deg_lst[src]), 1);
                __sync_fetch_and_add(&(deg_lst[dst]), 1);
            }
            // PrefixSum.
            {
                InclusivePrefixSumOMP(histogram, &off.front() + 1, num_vertices, [&deg_lst](uint32_t it) {
                    return deg_lst[it];
                }, max_omp_threads);
                size_t task_num = num_vertices + 1;
                size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
                auto it_beg = avg * tid;
                auto it_end = min(avg * (tid + 1), task_num);
                memcpy(cur_write_off + it_beg, &off.front() + it_beg, sizeof(uint32_t) * (it_end - it_beg));
            }
#pragma omp single
            {
                adj_lst = (int32_t *) malloc(sizeof(int32_t) * off[off.size() - 1]);
                log_info("before csr transform time: %.9lf s", convert_timer.elapsed());
            }
            // Scatter.
#pragma omp for
            for (uint32_t i = 0; i < num_edges; i++) {
                auto src = edge_lst[i].first;
                auto dst = edge_lst[i].second;

                auto old_offset = __sync_fetch_and_add(&(cur_write_off[src]), 1);
                adj_lst[old_offset] = dst;

                old_offset = __sync_fetch_and_add(&(cur_write_off[dst]), 1);
                adj_lst[old_offset] = src;
            }
#pragma omp single
            {
                log_info("before sort time: %.9lf s", convert_timer.elapsed());
            }
            // 5th: sort each local ranges
#pragma omp for schedule(dynamic)
            for (auto u = 0; u < deg_lst.size(); u++) {
                assert(cur_write_off[u] == off[u + 1]);
                sort(adj_lst + off[u], adj_lst + off[u + 1]);
            }
        }
        log_info("edge list to csr time: %.9lf s", convert_timer.elapsed());

        graph_t g;
        g.n = num_vertices;
        assert(off[off.size() - 1] == 2 * num_edges);
        g.m = off[off.size() - 1];
        g.adj = adj_lst;
        g.num_edges = &off.front();
        size_t tc_cnt = 0;
#pragma omp parallel for schedule(dynamic, 6000) reduction(+:tc_cnt)
        for (auto i = 0u; i < g.m; i++)
            ComputeSupport(&g, tc_cnt, i);
        log_info("graph :%lld, %lld", g.n, g.m);
        log_info("There are %zu triangles in the input graph.", tc_cnt / 3);
    }
}