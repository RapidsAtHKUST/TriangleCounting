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