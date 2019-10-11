#pragma once

#include <new>
#include <malloc.h>
#include "omp.h"

#include "../tc_utils.h"
#include "boolarray.h"

inline uint32_t get_log_size(int x) {
    int cnt = 0;
    for (; x > 0; cnt++) {
        x >>= 1;
    }
    return cnt;
}

inline uint32_t get_part_size(int i) {
    return i == 0 ? 0 : 1 << (get_log_size(i) - 1);
}

class RadixFilter {
    BoolArray<uint64_t> psum_occupied_bool_;
    graph_t *g_;
    uint32_t *row_ptrs_end;
    int radix_val_;
public:
    explicit RadixFilter(graph_t *g, uint32_t *row_ptrs_end) : g_(g), row_ptrs_end(row_ptrs_end), radix_val_(-1) {}

    void Construct(int u) {
        auto deg = row_ptrs_end[u + 1] - g_->row_ptrs[u];
//        if (deg > 0) {        // assume deg > 0
        // 1: Histogram.
        constexpr int heuristic_factor = 64;
        auto partition_size = get_part_size(deg) * heuristic_factor;
        radix_val_ = partition_size - 1;
        psum_occupied_bool_ = BoolArray<uint64_t>(partition_size);
        auto radix_val = partition_size - 1;
        for (auto j = g_->row_ptrs[u]; j < row_ptrs_end[u + 1]; j++) {
            auto v = g_->adj[j];
            auto bucket_id = v & radix_val;
            psum_occupied_bool_.set(bucket_id);
        }
//        }
    }

    bool PossibleExist(int v) {
        auto bucket_id = v & radix_val_;
        return psum_occupied_bool_.get(bucket_id);
    }
};


class RadixSet {
    vector<int> psum_arr_;
    BoolArray<uint64_t> psum_occupied_bool_;
    vector<int> tmp_;
    vector<int> hash_table_;
    graph_t *g_;
    uint32_t *row_ptrs_end;

public:
    explicit RadixSet(graph_t *g, uint32_t *row_ptrs_end) : g_(g), row_ptrs_end(row_ptrs_end) {
        psum_arr_.reserve(1024 * 1204 * 2);
        hash_table_.reserve(1024 * 1204 * 2);
    }

    void Construct(int u) {
        auto deg = row_ptrs_end[u + 1] - g_->row_ptrs[u];
        if (deg > 0) {
            // 1: Histogram.
            constexpr int heuristic_factor = 64;
            auto partition_size = get_part_size(deg) * heuristic_factor;
            psum_arr_.resize(partition_size + 1);
            memset(&psum_arr_.front(), 0, psum_arr_.size() * sizeof(int));
            hash_table_.resize(deg);
            psum_occupied_bool_ = BoolArray<uint64_t>(psum_arr_.size());
            auto radix_val = partition_size - 1;
            for (auto j = g_->row_ptrs[u]; j < row_ptrs_end[u + 1]; j++) {
                auto v = g_->adj[j];
                auto bucket_id = v & radix_val;
                psum_arr_[bucket_id + 1]++;
                psum_occupied_bool_.set(bucket_id);
            }
            // 2: PrefixSum.
            for (auto i = 0u; i < partition_size; i++) {
                psum_arr_[i + 1] += psum_arr_[i];
            }
            // 3: Scatter.
            tmp_.resize(psum_arr_.size());
            memcpy(&tmp_.front(), &psum_arr_.front(), sizeof(int) * psum_arr_.size());
            for (auto j = g_->row_ptrs[u]; j < row_ptrs_end[u + 1]; j++) {
                auto v = g_->adj[j];
                auto bucket_id = v & radix_val;
                hash_table_[tmp_[bucket_id]++] = v;
            }
        } else {
            psum_arr_.clear();
        }
    }

    bool Exist(int v) {
        if (psum_arr_.empty())return false;   // Assume not zero deg.

        auto partition_size = psum_arr_.size() - 1;

        auto radix_val = partition_size - 1;
        auto bucket_id = v & radix_val;

        if (psum_occupied_bool_.get(bucket_id)) {
            for (auto j = psum_arr_[bucket_id]; j < psum_arr_[bucket_id + 1]; j++) {
                if (hash_table_[j] == v) {
                    return true;
                }
            }
        }
        return false;
    }
};