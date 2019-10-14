#pragma once

#include "tc_utils.h"
#include "util/radix_hash_map.h"
#include "util/libpopcnt.h"
#include "util/set_inter_cnt_utils.h"

#define MAX_PACK_NUM (65536)

inline size_t CountTriBMPAndMergeWithPack(graph_t &g, int max_omp_threads) {
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

    int to_pack_num = min<int>(g.n, MAX_PACK_NUM);
    vector<vector<uint16_t>> word_indexes(to_pack_num);   // 65536 * bits-sizeof(word)
    vector<vector<word_t>> words(to_pack_num);

#pragma omp parallel num_threads(max_omp_threads)
    {
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
        for (auto u = 0u; u < to_pack_num; u++) {
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
            static thread_local BoolArray<word_t> bitmap(threshold);
            static thread_local vector<word_t> buffer(threshold / word_in_bits);
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
            auto du = row_ptrs_end[u + 1] - row_ptrs_beg[u];

            for (auto edge_idx = g.row_ptrs[u]; edge_idx < row_ptrs_end[u + 1]; edge_idx++) {
                auto v = g.adj[edge_idx];
                auto cn_count = 0;
                // First Range.
                if (v < to_pack_num) {
                    auto num_words_v = word_indexes[v].size();
                    for (size_t i = 0; i < num_words_v; i++) {
#ifdef WORKLOAD_STAT
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

                // Second Range.
                auto dv = row_ptrs_end[v + 1] - row_ptrs_beg[v];
                if (du > 0 && dv > 0) {
#ifdef __AVX2__
                    cn_count += SetInterCntAVX2Detail(&g, row_ptrs_beg[u], row_ptrs_end[u + 1], row_ptrs_beg[v],
                                                           row_ptrs_end[v + 1]);
#elif defined(__SSE4_1__)
                    cn_count += SetInterCntSSE4Detail(&g, row_ptrs_beg[u], row_ptrs_end[u + 1], row_ptrs_beg[v],
                                                      row_ptrs_end[v + 1]);
#else
                    cn_count += SetIntersectionScalarCntDetail(&g, row_ptrs_beg[u], row_ptrs_end[u + 1],
                                                               row_ptrs_beg[v], row_ptrs_end[v + 1]);
#endif
                }
                tc_cnt += cn_count;
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
    log_info("Workload (large-deg vid in [0, %d]): %s, avg: %s", threshold,
             FormatWithCommas(workload_large_deg).c_str(),
             FormatWithCommas(workload_large_deg / (g.m / 2)).c_str());
#endif
    return tc_cnt;
}

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
                                          min<int>(threshold, u)) - g.adj;
            max_d = max<int>(max_d, row_ptrs_end[u + 1] - g.row_ptrs[u]);
        }
#pragma omp single
        {
            log_info("finish init row_ptrs_end, max d: %d, time: %.9lfs", max_d, tc_timer.elapsed());
        }

        // Construct Words for Range [0, 32768), at most 512 words (given each word 64 bits)
#pragma omp for schedule(dynamic, 100)
        for (auto u = 0u; u < to_pack_num; u++) {
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
            static thread_local BoolArray<word_t> bitmap(threshold);
            if (u < to_pack_num) {
                for (size_t i = 0; i < word_indexes[u].size(); i++) {
//                    assert(word_indexes[u][i] < threshold / word_in_bits);
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
                    for (size_t i = 0; i < word_indexes[v].size(); i++) {
#ifdef WORKLOAD_STAT
                        workload_large_deg++;
#endif
                        word_t word = bitmap.getWord(word_indexes[v][i]) & words[v][i];
                        cn_count += popcnt(&word, sizeof(word_t));
                    }
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
#ifdef WORKLOAD_STAT
    log_info("Workload: %s, avg: %s", FormatWithCommas(workload).c_str(),
             FormatWithCommas(workload / (g.m / 2)).c_str());
    log_info("Workload (large-deg vid in [0, %d]): %s, avg: %s", threshold,
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