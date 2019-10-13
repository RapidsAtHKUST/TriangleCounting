#pragma once

#include <vector>
#include <cstdint>

#include <omp.h>

#include "local_buffer.h"

using namespace std;

#define CACHE_LINE_ENTRY (16)
#define LOCAL_BUDGET (8*1024*1024)

template<typename T>
void MemSetOMP(T *arr, int val, size_t size, size_t tid, size_t max_omp_threads) {
    size_t task_num = size;
    size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
    auto it_beg = avg * tid;
    auto it_end = min(avg * (tid + 1), task_num);
    memset(arr + it_beg, val, sizeof(T) * (it_end - it_beg));
#pragma omp barrier
}

template<typename T>
void MemCpyOMP(T *dst, T *src, size_t size, size_t tid, size_t max_omp_threads) {
    size_t task_num = size;
    size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
    auto it_beg = avg * tid;
    auto it_end = min(avg * (tid + 1), task_num);
    memcpy(dst + it_beg, src + it_beg, sizeof(T) * (it_end - it_beg));
#pragma omp barrier
}

/*
 * InclusivePrefixSumOMP: General Case Inclusive Prefix Sum
 * histogram: is for cache-aware thread-local histogram purpose
 * output: should be different from the variables captured in function object f
 * size: is the original size for the flagged prefix sum
 * f: requires it as the parameter, f(it) return the histogram value of that it
 */
template<typename T, typename F>
void InclusivePrefixSumOMP(vector<uint32_t> &histogram, T *output, size_t size, F f, int omp_num_threads) {
#pragma omp single
    {
        histogram = vector<uint32_t>((omp_num_threads + 1) * CACHE_LINE_ENTRY, 0);
    }
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
        for (auto local_tid = 0; local_tid < omp_num_threads; local_tid++) {
            auto local_histogram_idx = (local_tid + 1) * CACHE_LINE_ENTRY;
            auto prev_histogram_idx = (local_tid) * CACHE_LINE_ENTRY;
            histogram[local_histogram_idx] += histogram[prev_histogram_idx];
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
 * FlagPrefixSumOMP: special case of InclusivePrefixSumOMP
 */
template<typename T, typename F>
void FlagPrefixSumOMP(vector<uint32_t> &histogram, T *output, size_t size, F f, int omp_num_threads) {
    InclusivePrefixSumOMP(histogram, output, size, [&f](size_t it) {
        return f(it) ? 1 : 0;
    }, omp_num_threads);
}

/*
 * SelectNotFOMP: selection primitive
 * !f(it) returns selected
 */
template<typename T, typename OFF, typename F>
void SelectNotFOMP(vector<uint32_t> &histogram, T *output, T *input,
                   OFF *relative_off, size_t size, F f, int omp_num_threads) {
    FlagPrefixSumOMP(histogram, relative_off, size, f, omp_num_threads);
#pragma omp for
    for (size_t i = 0u; i < size; i++) {
        if (!(f(i))) {
            auto off = i - relative_off[i];
            output[off] = input[i];
        }
    }
}

/*
 * Require an output array,
 * f: is the property for the bucket ID, given an index on the input array
 * Inefficient when there are lots of contentions because of atomic operations
 */
template<typename T, typename OFF, typename F>
void BucketSort(vector<uint32_t> &histogram, T *&input, T *&output,
                OFF *&cur_write_off, OFF *&bucket_ptrs,
                size_t size, int32_t num_buckets, F f, int max_omp_threads, Timer *timer = nullptr) {
    int tid = omp_get_thread_num();

    // Populate.
#pragma omp single
    {
        bucket_ptrs = (OFF *) malloc(sizeof(OFF) * (num_buckets + 1));
        cur_write_off = (OFF *) malloc(sizeof(OFF) * (num_buckets + 1));
        cur_write_off[0] = 0;
    }
    MemSetOMP(bucket_ptrs, 0, num_buckets + 1, tid, max_omp_threads);

    // Histogram.
    auto local_buf = (uint8_t *) calloc(num_buckets, sizeof(uint8_t));
#pragma omp for
    for (size_t i = 0u; i < size; i++) {
        auto src = f(i);
        local_buf[src]++;
        if (local_buf[src] == 0xff) {
            __sync_fetch_and_add(&bucket_ptrs[src], 0xff);
            local_buf[src] = 0;
        }
    }
    for (size_t i = 0; i < num_buckets; i++) {
        if (local_buf[i] != 0) {
            __sync_fetch_and_add(&bucket_ptrs[i], local_buf[i]);
        }
    }
#pragma omp barrier
    free(local_buf);

    InclusivePrefixSumOMP(histogram, cur_write_off + 1, num_buckets, [&bucket_ptrs](uint32_t it) {
        return bucket_ptrs[it];
    }, max_omp_threads);
    MemCpyOMP(bucket_ptrs, cur_write_off, num_buckets + 1, tid, max_omp_threads);

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

template<typename T, typename OFF, typename F>
void BucketSortSmallBuckets(vector<uint32_t> &histogram, T *&input, T *&output,
                            OFF *&cur_write_off, OFF *&bucket_ptrs,
                            size_t size, int32_t num_buckets, F f, int max_omp_threads, Timer *timer = nullptr) {
    int tid = omp_get_thread_num();
    using BufT= LocalWriteBuffer<T, uint32_t>;
    auto cap = max<int>(CACHE_LINE_ENTRY, LOCAL_BUDGET / num_buckets / sizeof(T));
    auto bucket_write_buffers = (BufT *) malloc(num_buckets * sizeof(BufT));
    auto bucket_buffers = (T *) malloc(cap * num_buckets * sizeof(T));
    auto counter = (int *) calloc(num_buckets, sizeof(int));
    // Populate.
#pragma omp single
    {
        log_info("Mem Size Buckets: %zu, Bucket#: %d", cap * num_buckets * sizeof(T) * max_omp_threads, num_buckets);
        bucket_ptrs = (uint32_t *) malloc(sizeof(OFF) * (num_buckets + 1));
        cur_write_off = (uint32_t *) malloc(sizeof(OFF) * (num_buckets + 1));
        cur_write_off[0] = 0;
    }
    {
        size_t task_num = num_buckets + 1;
        size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
        auto it_beg = avg * tid;
        auto it_end = min(avg * (tid + 1), task_num);
        memset(bucket_ptrs + it_beg, 0, sizeof(OFF) * (it_end - it_beg));
#pragma omp barrier
    }
    // Histogram.
#pragma omp for
    for (size_t i = 0u; i < size; i++) {
        counter[f(i)]++;
    }
    for (auto i = 0; i < num_buckets; i++) {
        if (counter[i] > 0)
            __sync_fetch_and_add(&(bucket_ptrs[i]), counter[i]);
    }
#pragma omp barrier
    InclusivePrefixSumOMP(histogram, cur_write_off + 1, num_buckets, [&bucket_ptrs](uint32_t it) {
        return bucket_ptrs[it];
    }, max_omp_threads);

    {
        size_t task_num = num_buckets + 1;
        size_t avg = (task_num + max_omp_threads - 1) / max_omp_threads;
        auto it_beg = avg * tid;
        auto it_end = min(avg * (tid + 1), task_num);
        memcpy(bucket_ptrs + it_beg, cur_write_off + it_beg, sizeof(OFF) * (it_end - it_beg));
#pragma omp barrier
    }
#pragma omp single
    {
        if (timer != nullptr)log_info("Before Scatter, Time: %.9lfs", timer->elapsed());
    }
    for (auto i = 0; i < num_buckets; i++) {
        bucket_write_buffers[i] = BufT(bucket_buffers + cap * i, cap, output, &cur_write_off[i]);
    }
#pragma omp barrier
    // Scatter.
#pragma omp for
    for (size_t i = 0u; i < size; i++) {
        auto element = input[i];
        auto bucket_id = f(i);
        bucket_write_buffers[bucket_id].push(element);
    }
    for (auto i = 0; i < num_buckets; i++) {
        bucket_write_buffers[i].submit_if_possible();
    }
#pragma omp barrier
#pragma omp single
    {
        if (timer != nullptr)log_info("Before Sort, Time: %.9lfs", timer->elapsed());
    }

    free(counter);
    free(bucket_buffers);
    free(bucket_write_buffers);
#pragma omp barrier
}