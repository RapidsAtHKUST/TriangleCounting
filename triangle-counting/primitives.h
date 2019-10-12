#pragma once

#include <vector>
#include <cstdint>

#include <omp.h>

using namespace std;

#define CACHE_LINE_ENTRY (16)

/*
 * Both Primitives Are Inclusive
 * histogram is for cache-aware thread-local histogram purpose
 * output should be different from the variables captured in function object f
 * size is the original size for the flagged prefix sum
 */

template<typename F>
void FlagPrefixSumOMP(vector<uint32_t> &histogram, uint32_t *output, size_t size, F f, int omp_num_threads) {
    static thread_local int tid = omp_get_thread_num();
    // 1st Pass: Histogram.
    auto avg = size / omp_num_threads;
    auto it_beg = avg * tid;
    auto histogram_idx = (tid + 1) * CACHE_LINE_ENTRY;
    histogram[histogram_idx] = 0;
    auto it_end = tid == omp_num_threads - 1 ? size : avg * (tid + 1);
    auto prev = 0u;
    for (auto it = it_beg; it < it_end; it++) {
        if (f(it)) {
            histogram[histogram_idx]++;
            output[it] = prev + 1;
        } else {
            output[it] = prev;
        }
        prev = output[it];
    }
#pragma omp barrier

    // 2nd Pass: single-prefix-sum & Add previous sum.
#pragma omp single
    {
        for (auto tid = 0; tid < omp_num_threads; tid++) {
            auto histogram_idx = (tid + 1) * CACHE_LINE_ENTRY;
            auto prev_histogram_idx = (tid) * CACHE_LINE_ENTRY;
            histogram[histogram_idx] += histogram[prev_histogram_idx];
        }
    }
    {
        auto prev_sum = histogram[tid * CACHE_LINE_ENTRY];
        for (auto it = it_beg; it < it_end; it++) {
            output[it] += prev_sum;
        }
#pragma omp barrier
    }
}

template<typename F>
void InclusivePrefixSumOMP(vector<uint32_t> &histogram, uint32_t *output, size_t size, F f, int omp_num_threads) {
    static thread_local int tid = omp_get_thread_num();
    // 1st Pass: Histogram.
    auto avg = size / omp_num_threads;
    auto it_beg = avg * tid;
    auto histogram_idx = (tid + 1) * CACHE_LINE_ENTRY;
    histogram[histogram_idx] = 0;
    auto it_end = tid == omp_num_threads - 1 ? size : avg * (tid + 1);
    auto prev = 0u;
    for (auto it = it_beg; it < it_end; it++) {
        auto value = f(it);
        histogram[histogram_idx] += value;
        prev += value;
        output[it] = prev;
    }
#pragma omp barrier

    // 2nd Pass: single-prefix-sum & Add previous sum.
#pragma omp single
    {
        for (auto tid = 0; tid < omp_num_threads; tid++) {
            auto histogram_idx = (tid + 1) * CACHE_LINE_ENTRY;
            auto prev_histogram_idx = (tid) * CACHE_LINE_ENTRY;
            histogram[histogram_idx] += histogram[prev_histogram_idx];
        }
    }
    {
        auto prev_sum = histogram[tid * CACHE_LINE_ENTRY];
        for (auto it = it_beg; it < it_end; it++) {
            output[it] += prev_sum;
        }
#pragma omp barrier
    }
}

/*
 * Require an output array,
 * f: is the property for the bucket ID, given an index on the input array
 * Inefficient when there are lots of contentions because of atomic operations
 */
template<typename T, typename O, typename F>
void BucketSort(vector<uint32_t> &histogram, T *&input, T *&output,
                O *&cur_write_off, O *&bucket_ptrs,
                size_t size, int32_t num_buckets, F f, int max_omp_threads, Timer *timer = nullptr) {
    int tid = omp_get_thread_num();

    // Populate.
#pragma omp single
    {
        bucket_ptrs = (uint32_t *) malloc(sizeof(O) * (num_buckets + 1));
        cur_write_off = (uint32_t *) malloc(sizeof(O) * (num_buckets + 1));
        cur_write_off[0] = 0;
    }
    {
        size_t task_num = num_buckets + 1;
        size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
        auto it_beg = avg * tid;
        auto it_end = min(avg * (tid + 1), task_num);
        memset(bucket_ptrs + it_beg, 0, sizeof(O) * (it_end - it_beg));
#pragma omp barrier
    }
    // Histogram.
#pragma omp for
    for (size_t i = 0u; i < size; i++) {
        __sync_fetch_and_add(&(bucket_ptrs[f(i)]), 1);
    }
    InclusivePrefixSumOMP(histogram, cur_write_off + 1, num_buckets, [&bucket_ptrs](uint32_t it) {
        return bucket_ptrs[it];
    }, max_omp_threads);

    {
        size_t task_num = num_buckets + 1;
        size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
        auto it_beg = avg * tid;
        auto it_end = min(avg * (tid + 1), task_num);
        memcpy(bucket_ptrs + it_beg, cur_write_off + it_beg, sizeof(uint32_t) * (it_end - it_beg));
#pragma omp barrier
    }
#pragma omp single
    {
        if (timer != nullptr)log_info("Before Scatter, Time: %.9lfs", timer->elapsed());
    }
    // Scatter.
#pragma omp for
    for (size_t i = 0u; i < size; i++) {
        auto element = input[i];
        auto bucket_id = f(i);
        auto old_offset = __sync_fetch_and_add(&(cur_write_off[bucket_id]), 1);
        output[old_offset] = element;
    }
#pragma omp single
    {
        if (timer != nullptr)log_info("Before Sort, Time: %.9lfs", timer->elapsed());
    }
#pragma omp barrier
}