```cpp
#ifdef VERIFY_TABLE
        int32_t tmp = -1;
        for (uint32_t i = 0; i < num_edges; i++) {
            tmp = max(tmp, max(edge_lst[i].first, edge_lst[i].second));
        }
        vector<unordered_set<int32_t >> tables(tmp + 1);
        for (auto i = 0u; i < num_edges; i++) {
            auto edge = edge_lst[i];
            if (edge.first != edge.second) {
                tables[min(edge.first, edge.second)].emplace(max(edge.first, edge.second));
            }
        }
        size_t gt_size = 0;
        for (auto &table:tables) {
            gt_size += table.size();
        }
        log_info("GT edge#: %zu", gt_size);
#endif
```

```cpp
template<typename T, typename I>
T RemoveDuplicates(pair<T, T> *&edge_lst, I &num_edges) {
    using Edge = pair<T, T>;
    Timer timer;
    T max_node_id = 0;
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
//#define BASELINE_SORT
#ifndef BASELINE_SORT
    parasort(num_edges, edge_lst, [](const Edge &left, const Edge &right) {
        if (left.first == right.first) {
            return left.second < right.second;
        }
        return left.first < right.first;
    }, max_omp_threads);
#else
    sort(edge_lst, edge_lst + num_edges,
         [](const Edge &left, const Edge &right) {
             if (left.first == right.first) {
                 return left.second < right.second;
             }
             return left.first < right.first;
         });
#endif
    log_info("Finish Sort, %.9lfs", timer.elapsed());

//#define NAIVE_REMOVE_DUPLICATE
#ifndef NAIVE_REMOVE_DUPLICATE
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
#else
    vector<Edge> edges2;

    auto last_u = -1;
    auto last_v = -1;
    for (auto i = 0u; i < num_edges; i++) {
        auto edge = edge_lst[i];
        if (edge.first != last_u) {
            last_u = edge.first;
            last_v = edge.second;
            if (edge.first != edge.second)
                edges2.emplace_back(edge);
        } else {
            if (edge.second != last_v) {
                last_v = edge.second;
                if (edge.first != edge.second)
                    edges2.emplace_back(edge);
            }
        }
    }
    edge_lst = &edges2.front();
    num_edges = edges2.size();
#endif
    log_info("New # of edges: %zu, Elapsed: %.9lfs", num_edges, timer.elapsed());
    return max_node_id;
}
```