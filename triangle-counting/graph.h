#pragma once

#include <cstdint>

struct graph_t {
    long n;
    long m;

    int32_t *adj;
    uint32_t *row_ptrs;
};



