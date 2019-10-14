```cpp

inline size_t CountTriBMPWithPack(graph_t &g, int max_omp_threads) {
    Timer tc_timer;
    int max_d = 0;
    size_t tc_cnt = 0;
    auto *row_ptrs_end = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
    auto *row_ptrs_beg = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
    size_t workload = 0;
    size_t workload_large_deg = 0;

    using word_t = uint64_t;
    constexpr int word_in_bits = sizeof(word_t) * 8;

    int to_pack_num = min<int>(g.n, MAX_PACK_NUM);
    vector<vector<uint16_t>> word_indexes(to_pack_num);   // 65536 * bits-sizeof(word)
    vector<vector<word_t>> words(to_pack_num);

#pragma omp parallel num_threads(max_omp_threads)
    {
        auto bits_vec = vector<bool>(g.n, false);
#pragma omp for reduction(max: max_d)
        for (auto u = 0u; u < g.n; u++) {
            row_ptrs_end[u + 1] = static_cast<uint32_t>(
                    lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1], u) - g.adj);
            row_ptrs_beg[u] = lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1],
                                          min<int>(FIRST_RANGE_SIZE, u)) - g.adj;
            max_d = max<int>(max_d, row_ptrs_end[u + 1] - g.row_ptrs[u]);
        }
#pragma omp single
        {
            log_info("finish init row_ptrs_end, max d: %d, time: %.9lfs", max_d, tc_timer.elapsed());
        }

        PackWords(g, row_ptrs_beg, to_pack_num, word_indexes, words, tc_timer);

#pragma omp for schedule(dynamic, 100) reduction(+:tc_cnt) reduction(+:workload)  reduction(+:workload_large_deg)
        for (auto u = 0u; u < g.n; u++) {
            // Set.
            static thread_local BoolArray<word_t> bitmap(FIRST_RANGE_SIZE);
            static thread_local vector<word_t> buffer(FIRST_RANGE_SIZE / word_in_bits);

            if (u < to_pack_num) {
                for (size_t i = 0; i < word_indexes[u].size(); i++) {
                    bitmap.setWord(word_indexes[u][i], words[u][i]);
                }
            } else {
                for (auto off = g.row_ptrs[u]; off < row_ptrs_beg[u]; off++) {
                    auto w = g.adj[off];
                    bitmap.set(w);
                }
            }

            for (auto edge_idx = row_ptrs_beg[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                bits_vec[v] = true;
            }

            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                auto cn_count = 0;
                if (v < to_pack_num) {
                    auto num_words_v = word_indexes[v].size();
                    for (size_t i = 0; i < num_words_v; i++) {
#ifdef WORKLOAD_STAT
                        workload++;
                        workload_large_deg++;
#endif
                        buffer[i] = bitmap.getWord(word_indexes[v][i]);
                    }
                    for (size_t i = 0; i < num_words_v; i++) {
                        buffer[i] &= words[v][i];
                    }
                    cn_count += popcnt(&buffer.front(), sizeof(word_t) * num_words_v);
                } else {
                    for (auto off = g.row_ptrs[v]; off < row_ptrs_beg[v]; off++) {
                        auto w = g.adj[off];
                        if (bitmap[w])cn_count++;
                    }
                }

                for (auto offset = row_ptrs_beg[v]; offset < row_ptrs_end[v + 1]; offset++) {
#ifdef WORKLOAD_STAT
                    workload++;
#endif
                    auto w = g.adj[offset];
                    if (bits_vec[w]) {
                        cn_count++;
                    }
                }
                tc_cnt += cn_count;
            }

            // Clear.
            for (auto edge_idx = row_ptrs_beg[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                bits_vec[v] = false;
            }
            bitmap.reset();
        }
    }
    free(row_ptrs_beg);
    free(row_ptrs_end);
    log_info("Forward cost: %.3lf s, Mem Usage: %d KB", tc_timer.elapsed(), getValue());
    log_info("Triangle Cnt: %'zu", tc_cnt);
#ifdef WORKLOAD_STAT
    log_info("Workload: %s, avg: %s", FormatWithCommas(workload).c_str(),
             FormatWithCommas(workload / (g.m / 2)).c_str());
    log_info("Workload (large-deg vid in [0, %d]): %s, avg: %s", FIRST_RANGE_SIZE,
             FormatWithCommas(workload_large_deg).c_str(),
             FormatWithCommas(workload_large_deg / (g.m / 2)).c_str());
#endif
    return tc_cnt;
}

inline size_t CountTriBMP(graph_t &g, int max_omp_threads) {
    Timer tc_timer;
    int max_d = 0;
    size_t tc_cnt = 0;
    auto *row_ptrs_end = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
    size_t workload = 0;
    size_t workload_large_deg = 0;
#pragma omp parallel num_threads(max_omp_threads)
    {
        // Bit vectors.
        auto bits_vec = vector<bool>(g.n, false);
#pragma omp for reduction(max: max_d)
        for (auto u = 0u; u < g.n; u++) {
            row_ptrs_end[u + 1] = static_cast<uint32_t>(
                    lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1], u) - g.adj);
            max_d = max<int>(max_d, row_ptrs_end[u + 1] - g.row_ptrs[u]);
        }
#pragma omp single
        {
            log_info("finish init row_ptrs_end, max d: %d", max_d);
            {
                stringstream ss;
                ss << pretty_print_array(g.adj + g.row_ptrs[0], row_ptrs_end[1] - g.row_ptrs[0]);
                ss << pretty_print_array(g.adj + g.row_ptrs[1000], row_ptrs_end[1001] - g.row_ptrs[1000]);
                log_info("Example v0, v1000: %s", ss.str().c_str());
            }
        };

#pragma omp for schedule(dynamic, 100) reduction(+:tc_cnt) reduction(+:workload)  reduction(+:workload_large_deg)
        for (auto u = 0u; u < g.n; u++) {
            // Set.
            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                bits_vec[v] = true;
            }
            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                auto cn_count = 0;
                for (auto offset = g.row_ptrs[v]; offset < row_ptrs_end[v + 1]; offset++) {
                    workload++;
                    auto w = g.adj[offset];
                    workload_large_deg += w < FIRST_RANGE_SIZE ? 1 : 0;

                    if (bits_vec[w]) {
                        cn_count++;
                    }
                }
                tc_cnt += cn_count;
            }
            // Clear.
            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                bits_vec[v] = false;
            }
        }
    }
    free(row_ptrs_end);
    log_info("Forward cost: %.3lf s, Mem Usage: %d KB",
             tc_timer.elapsed(), getValue());
    log_info("Triangle Cnt: %'zu", tc_cnt);
#ifdef WORKLOAD_STAT
    log_info("Workload: %s, avg: %s", FormatWithCommas(workload).c_str(),
             FormatWithCommas(workload / (g.m / 2)).c_str());
    log_info("Workload (large-deg vid in [0, %d]): %s, avg: %s", FIRST_RANGE_SIZE,
             FormatWithCommas(workload_large_deg).c_str(),
             FormatWithCommas(workload_large_deg / (g.m / 2)).c_str());
#endif
    return tc_cnt;
}

inline size_t CountTriRadixPartition(graph_t &g, int max_omp_threads) {
    Timer tc_timer;
    int max_d = 0;
    size_t tc_cnt = 0;
    auto *row_ptrs_end = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
#pragma omp parallel num_threads(max_omp_threads)
    {
        // Bit vectors.
        auto bits_vec = vector<bool>(g.n, false);
#pragma omp for reduction(max: max_d)
        for (auto u = 0u; u < g.n; u++) {
            row_ptrs_end[u + 1] = static_cast<uint32_t>(
                    lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1], u) - g.adj);
            max_d = max<int>(max_d, row_ptrs_end[u + 1] - g.row_ptrs[u]);
        }
#pragma omp single
        log_info("finish init row_ptrs_end, max d: %d", max_d);

#pragma omp for schedule(dynamic, 100) reduction(+:tc_cnt)
        for (auto u = 0u; u < g.n; u++) {
            // Set.
            static thread_local RadixSet radix_set(&g, row_ptrs_end);
            radix_set.Construct(u);
            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                auto cn_count = 0;
                for (auto offset = g.row_ptrs[v]; offset < row_ptrs_end[v + 1]; offset++) {
                    auto w = g.adj[offset];
                    if (radix_set.Exist(w)) {
                        cn_count++;
                    }
                }
                tc_cnt += cn_count;
            }
        }
    }
    free(row_ptrs_end);
    log_info("Forward cost: %.3lf s, Mem Usage: %d KB",
             tc_timer.elapsed(), getValue());
    log_info("Triangle Cnt: %zu", tc_cnt);
    return tc_cnt;
}

```