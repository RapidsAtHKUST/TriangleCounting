#pragma once

#include "tc_utils.h"
#include "util/radix_hash_map.h"
#include "util/libpopcnt.h"

inline size_t CountTriBMPWithPack(graph_t &g, int max_omp_threads) {
    Timer tc_timer;
    int max_d = 0;
    size_t tc_cnt = 0;
    auto *row_ptrs_end = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
    auto *row_ptrs_beg = (uint32_t *) malloc(sizeof(uint32_t) * (g.n + 1));
    size_t workload = 0;
    size_t workload_large_deg = 0;
    int threshold = 32768; // word-packing for [0, 32768), otherwise normal table (using bitmap)

    using word_t = uint64_t;
    constexpr int word_in_bits = sizeof(word_t) * 8;
    vector<vector<uint16_t>> word_indexes(g.n);   // 65536 * bits-sizeof(word)
    vector<vector<word_t>> words(g.n);

#pragma omp parallel num_threads(max_omp_threads)
    {
        auto bits_vec = vector<bool>(g.n, false);
#pragma omp for reduction(max: max_d)
        for (auto u = 0u; u < g.n; u++) {
            row_ptrs_end[u + 1] = static_cast<uint32_t>(
                    lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1], u) - g.adj);
            row_ptrs_beg[u] = lower_bound(g.adj + g.row_ptrs[u], g.adj + g.row_ptrs[u + 1],
                                          min<int>(threshold, u)) - g.adj;
            max_d = max<int>(max_d, row_ptrs_end[u + 1] - g.row_ptrs[u]);
        }
#pragma omp single
        {
            log_info("finish init row_ptrs_end, max d: %d, time: %.9lfs", max_d, tc_timer.elapsed());
        }

        // Construct Words for Range [0, 32768), at most 512 words (given each word 64 bits)
#pragma omp for schedule(dynamic, 100)
        for (auto u = 0u; u < g.n; u++) {
            auto prev_blk_id = -1;
            for (auto off = g.row_ptrs[u]; off < row_ptrs_beg[u]; off++) {
                auto v = g.adj[off];
                int cur_blk_id = v / word_in_bits;
                if (cur_blk_id != prev_blk_id) {
                    prev_blk_id = cur_blk_id;
                    word_indexes[u].emplace_back(cur_blk_id);
                    words[u].emplace_back(0);
                }
                words[u].back() |= static_cast<word_t>(1u) << (v % word_in_bits);
            }
        }
#pragma omp single
        {
            log_info("Finish Indexing: %.9lfs", tc_timer.elapsed());
        }

#pragma omp for schedule(dynamic, 100) reduction(+:tc_cnt) reduction(+:workload)  reduction(+:workload_large_deg)
        for (auto u = 0u; u < g.n; u++) {
            // Set.
            BoolArray<word_t> bitmap(threshold);
            BoolArray<uint32_t> word_existence(threshold / word_in_bits);
            for (size_t i = 0; i < word_indexes[u].size(); i++) {
                bitmap.setWord(word_indexes[u][i], words[u][i]);
                word_existence.set(word_indexes[u][i]);
            }
            for (auto edge_idx = row_ptrs_beg[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                bits_vec[v] = true;
            }
            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                auto cn_count = 0;
                for (size_t i = 0; i < word_indexes[v].size(); i++) {
#ifdef WORKLOAD_STAT
                    workload_large_deg++;
#endif
                    if (word_existence.get(word_indexes[v][i])) {
                        word_t word = bitmap.getWord(word_indexes[v][i]) & words[v][i];
                        cn_count += popcnt(&word, sizeof(word_t));
                    }
                }
                for (auto offset = row_ptrs_beg[v]; offset < row_ptrs_end[v + 1]; offset++) {
#ifdef WORKLOAD_STAT
                    workload++;
#endif
                    auto w = g.adj[offset];
#ifdef WORKLOAD_STAT
                    workload_large_deg += w < threshold ? 1 : 0;
#endif
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
            for (size_t i = 0; i < word_indexes[u].size(); i++) {
                bitmap.setWord(word_indexes[u][i], 0);
            }
            word_existence.reset();
        }
    }
    free(row_ptrs_beg);
    free(row_ptrs_end);
    log_info("Forward cost: %.3lf s, Mem Usage: %d KB", tc_timer.elapsed(), getValue());
    log_info("Triangle Cnt: %'zu", tc_cnt);
#ifdef WORKLOAD_STAT
    log_info("Workload: %s, avg: %s", FormatWithCommas(workload).c_str(),
             FormatWithCommas(workload / (g.m / 2)).c_str());
    log_info("Workload (large-deg vid in [0, %d]): %s, avg: %s", threshold,
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
    int threshold = 32768;
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
                    workload_large_deg += w < threshold ? 1 : 0;

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
    log_info("Workload: %s, avg: %s", FormatWithCommas(workload).c_str(),
             FormatWithCommas(workload / (g.m / 2)).c_str());
    log_info("Workload (large-deg vid in [0, %d]): %s, avg: %s", threshold,
             FormatWithCommas(workload_large_deg).c_str(),
             FormatWithCommas(workload_large_deg / (g.m / 2)).c_str());
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