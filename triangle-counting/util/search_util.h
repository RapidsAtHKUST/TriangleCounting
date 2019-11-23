#pragma once

#include <x86intrin.h>

#include <cstdint>

template<typename T, typename OFF>
uint32_t LinearSearch(T *array, OFF offset_beg, OFF offset_end, T val) {
    // linear search fallback
    for (auto offset = offset_beg; offset < offset_end; offset++) {
        if (array[offset] >= val) {
            return offset;
        }
    }
    return offset_end;
}

template<typename T, typename OFF>
uint32_t BranchFreeBinarySearch(T *a, OFF offset_beg, OFF offset_end, T x) {
    OFF n = offset_end - offset_beg;
    const T *base = a + offset_beg;
    while (n > 1) {
        auto half = n / 2;
        __builtin_prefetch(base + half / 2, 0, 0);
        __builtin_prefetch(base + half + half / 2, 0, 0);
        base = (base[half] < x) ? base + half : base;
        n -= half;
    }
    return (*base < x) + base - a;
}

// Assuming (offset_beg != offset_end)
template<typename T, typename OFF>
uint32_t GallopingSearch(T *array, OFF offset_beg, OFF offset_end, T val) {
    if (array[offset_end - 1] < val) {
        return offset_end;
    }
    // galloping
    if (array[offset_beg] >= val) {
        return offset_beg;
    }
    if (array[offset_beg + 1] >= val) {
        return offset_beg + 1;
    }
    if (array[offset_beg + 2] >= val) {
        return offset_beg + 2;
    }

    OFF jump_idx = 4u;
    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BranchFreeBinarySearch(array, (jump_idx >> 1u) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1u;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BranchFreeBinarySearch(array, (jump_idx >> 1u) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}

#include <x86intrin.h>

#ifdef __AVX2__

template<typename OFF>
OFF LinearSearchAVX2(const int32_t *array, OFF offset_beg, OFF offset_end, int32_t val) {
    // linear search fallback
    __m256i pivot_element = _mm256_set1_epi32(val);
    for (; offset_beg + 8 < offset_end; offset_beg += 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        uint32_t mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) { return offset_beg + (_popcnt32(mask) >> 2); }
    }

    // offset_end - offset_beg <= 8
    auto left_size = offset_end - offset_beg;
    __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
    __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
    int mask = _mm256_movemask_epi8(cmp_res);
    int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
    mask &= cmp_mask;
    if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    return offset_end;
}

template<typename OFF>
OFF BinarySearchForGallopingSearchAVX2(const int32_t *array, OFF offset_beg, OFF offset_end, int32_t val) {
    while (offset_end - offset_beg >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_beg) + offset_end) / 2);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(mid + 1) + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(offset_beg) + mid) / 2], _MM_HINT_T0);
        if (array[mid] == val) {
            return mid;
        } else if (array[mid] < val) {
            offset_beg = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback, be careful with operator>> and operation+ priority
    __m256i pivot_element = _mm256_set1_epi32(val);
    for (; offset_beg + 8 < offset_end; offset_beg += 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        uint32_t mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) { return offset_beg + (_popcnt32(mask) >> 2); }
    }

    // offset_end - offset_beg <= 8
    auto left_size = offset_end - offset_beg;
    __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
    __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
    int mask = _mm256_movemask_epi8(cmp_res);
    int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
    mask &= cmp_mask;
    if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    return offset_end;
}

// Assuming size > 0
template<typename OFF>
uint32_t GallopingSearchAVX2(int *array, OFF offset_beg, OFF offset_end, int val) {
// Not necessary because of the linear search.
//    if (array[offset_end - 1] < val) {
//        return offset_end;
//    }

    // linear search
    __m256i pivot_element = _mm256_set1_epi32(val);
    if (offset_end - offset_beg >= 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        uint32_t mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) { return offset_beg + (_popcnt32(mask) >> 2); }
    } else {
        auto left_size = offset_end - offset_beg;
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    }

    // galloping, should add pre-fetch later
    auto jump_idx = 8u;
    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BranchFreeBinarySearch(array, (jump_idx >> 1) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BranchFreeBinarySearch(array, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}

#endif

#ifdef __AVX512F__

template<typename OFF>
uint32_t LinearSearchAVX512(int *array, OFF offset_beg, OFF offset_end, int val) {
    constexpr int parallelism = 16;
    __m512i pivot_element = _mm512_set1_epi32(val);
    for (; offset_beg + 15 < offset_end; offset_beg += parallelism) {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        if (mask != 0xffff) { return offset_beg + _mm_popcnt_u32(mask); }
    }
    if (offset_beg < offset_end) {
        auto left_size = offset_end - offset_beg;
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        __mmask16 cmp_mask = ((__mmask16) 0xffff) >> (16 - left_size);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + _mm_popcnt_u32(mask); }
    }
    return offset_end;
}

template<typename OFF>
uint32_t BinarySearchForGallopingSearchAVX512(const int *array, OFF offset_beg, OFF offset_end, int val) {
    while (offset_end - offset_beg >= 32) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_beg) + offset_end) / 2);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(mid + 1) + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(offset_beg) + mid) / 2], _MM_HINT_T0);
        if (array[mid] == val) {
            return mid;
        } else if (array[mid] < val) {
            offset_beg = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    constexpr int parallelism = 16;
    __m512i pivot_element = _mm512_set1_epi32(val);
    for (; offset_beg + 15 < offset_end; offset_beg += parallelism) {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        if (mask != 0xffff) { return offset_beg + _mm_popcnt_u32(mask); }
    }
    if (offset_beg < offset_end) {
        auto left_size = offset_end - offset_beg;
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        __mmask16 cmp_mask = ((__mmask16) 0xffff) >> (16 - left_size);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + _mm_popcnt_u32(mask); }
    }
    return offset_end;
}

template<typename OFF>
uint32_t GallopingSearchAVX512(int *array, OFF offset_beg, OFF offset_end, int val) {
    if (array[offset_end - 1] < val) {
        return offset_end;
    }

    // front peeking
    __m512i pivot_element = _mm512_set1_epi32(val);
    auto left_size = offset_end - offset_beg;
    if (left_size >= 16) {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        if (mask != 0xffff) { return offset_beg + _mm_popcnt_u32(mask); }
    } else {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        __mmask16 cmp_mask = ((__mmask16) 0xffff) >> (16 - left_size);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + _mm_popcnt_u32(mask); }
    }

    // galloping
    auto jump_idx = 16u;
    // pre-fetch
//    auto jump_times = 32 - _lzcnt_u32((offset_end - offset_beg) >> 4);
//    __m512i prefetch_idx = _mm512_set_epi64(16, 32, 64, 128, 256, 512, 1024, 2048);
//    __mmask8 mask = jump_times >= 8 ? (__mmask8) 0xff : (__mmask8) 0xff << (8 - jump_times);
//    _mm512_mask_prefetch_i64gather_ps(prefetch_idx, mask, array + offset_beg, 1, _MM_HINT_T0);

    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BinarySearchForGallopingSearchAVX512(array, (jump_idx >> 1) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BinarySearchForGallopingSearchAVX512(array, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}

#endif