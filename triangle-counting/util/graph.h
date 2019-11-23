#pragma once

#include <cstdint>

using row_ptr_t = uint64_t;

struct graph_t {
    long n;
    long m;

    int32_t *adj;
    row_ptr_t *row_ptrs;
};



