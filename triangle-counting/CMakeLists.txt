cmake_minimum_required(VERSION 3.5)
project(tri_cnt)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_compile_options(-Wall -O3 -g -march=native)
add_definitions(-DUSE_LOG)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else ()
    # icpc
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
        message(icpc)
        set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message(g++)
        set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif ()
endif ()

include_directories(.)
link_libraries(atomic)
add_executable(test_inplace_sort playground/playground_inplace_sort.cpp util/log.cpp)
add_executable(remove_duplicate playground/output_remove_duplicate.cpp util/log.cpp)
add_executable(tc triangle_counting_main.cpp util/log.cpp)

### Legacy Experimental Executables
### CSR-DODG-version
add_executable(tc-legacy-dodg triangle_counting_main_legacy.cpp util/log.cpp)
target_compile_definitions(tc-legacy-dodg PRIVATE WORKLOAD_STAT=1)
target_compile_definitions(tc-legacy-dodg PRIVATE DODG=1)

### CSR-normal-version
add_executable(tc-legacy triangle_counting_main_legacy.cpp util/log.cpp)
target_compile_definitions(tc-legacy PRIVATE WORKLOAD_STAT=1)
###