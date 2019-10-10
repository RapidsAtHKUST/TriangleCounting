add_compile_options(-Wall -O3 -g -march=native)
#add_executable(test_tbb test_tbb.cpp)
add_executable(test_reorder test_reorder.cpp ../reordering_utils/pkc_reorder.cpp ../reordering_utils/reorder_utils.cpp
        ../util/log.cpp ../util/graph.cpp ../util/md5.cpp ../util/file_system/file_util.h)

add_executable(test_bit_op test_bit_op ../util/log.cpp)

add_executable(test_inheritence test_inheritence.cpp ../util/log.cpp)

add_executable(clear_page_cache clear_page_cache.cpp)

add_executable(test_input test_input.cpp)