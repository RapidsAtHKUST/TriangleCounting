cmake_minimum_required(VERSION 3.6)
project(truss_decomp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH}) # To allow CMake to locate our Find*.cmake files
set(CMAKE_VERBOSE_MAKEFILE ON)

add_compile_options(-Wall -O3 -g -march=native)
add_definitions(-DUSE_LOG -DTIME_RESULTS)
option(USE_JEMALLOC "use jemalloc" ON)
option(PLAYGROUND "compile playground" ON)
option(USE_TCMALLOC "use tcmalloc" OFF)
option(BUILD_SERIAL "build serial codes" OFF)
option(LEGACY_PKT_ALL "build pkt all legacy codes" OFF)

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

#find_package(TBB REQUIRED)
#include_directories(${TBB_INCLUDE_DIRS})
#link_libraries(${TBB_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR}/../dependencies/libpopcnt)
include_directories(${CMAKE_SOURCE_DIR})

if (USE_TCMALLOC)
    # 1st: tcmalloc
    message("USE_TCMALLOC")
    find_package(Tcmalloc REQUIRED)
    if (Tcmalloc_FOUND)
        message("TCMALLOC FOUND: ${Tcmalloc_LIBRARY}")
        link_libraries(${Tcmalloc_LIBRARY})
    else ()
        message("NO TCMALLOC ")
    endif ()
elseif (USE_JEMALLOC)
    # 2nd: jemalloc
    message("USE_JEMALLOC")

    find_package(Jemalloc REQUIRED)
    if (JEMALLOC_FOUND)
        message("Link: ${JEMALLOC_LIBRARY}")
        link_libraries(${JEMALLOC_LIBRARY})
    else ()
        message("NO JEMALLOC ")
    endif ()
else ()
    message("Default Malloc")
    #    link_libraries("/usr/local/lib/libmimalloc.so")
endif ()

set(REORDER_UTILS reordering_utils/reorder_utils.cpp reordering_utils/pkc_reorder.cpp)
set(COMMON_UTILS util/graph.cpp util/stat.cpp util/log.cpp util/md5.cpp util/search_util.cpp)
add_library(common-utils ${REORDER_UTILS} ${COMMON_UTILS} extern_variables.cpp)
link_libraries(common-utils)
add_library(graph-pre-processing util/graph_preprocessing.cpp)

function(build_eid_ours name macros files)
    add_executable(${name} ${files} main.cpp iter_helper.cpp pkt_tri_cnt_intermediate.cpp util/graph_preprocessing.cpp)
    target_compile_definitions(${name} PRIVATE ${macros})
endfunction()

function(build_eval_tc name macros files)
    add_executable(${name} ${files} exp_varying_parameters/main_eval_tc.cpp iter_helper.cpp pkt_tri_cnt_intermediate.cpp)
    target_compile_definitions(${name} PRIVATE ${macros})
    target_link_libraries(${name} PRIVATE graph-pre-processing)
endfunction()

function(build_truss_algorithm name macros files)
    add_executable(${name} ${files} main.cpp iter_helper.cpp pkt_tri_cnt_intermediate.cpp)
    target_link_libraries(${name} PRIVATE graph-pre-processing)

    target_compile_definitions(${name} PRIVATE ${macros})
    target_compile_definitions(${name} PRIVATE "BMP_PROCESSED=1")
    target_compile_definitions(${name} PRIVATE "BMP_QUEUE=1")
    #    target_compile_definitions(${name} PRIVATE "LOCAL_STAT=1") # switch off when do not collect local stat
    #    target_compile_definitions(${name} PRIVATE "SUP_STAT=1") # switch off when do not collect local stat
endfunction()

function(build_varying_parameters name macros files)
    add_executable(${name} ${files} exp_varying_parameters/main_varying_gt_ept_rs.cpp
            iter_helper.cpp pkt_tri_cnt_intermediate.cpp pkt_merge.cpp)
    target_link_libraries(${name} PRIVATE graph-pre-processing)

    target_compile_definitions(${name} PRIVATE ${macros})
    target_compile_definitions(${name} PRIVATE "BMP_PROCESSED=1")
    target_compile_definitions(${name} PRIVATE "BMP_QUEUE=1")
    target_compile_definitions(${name} PRIVATE "PKT_INTER=1")
    #    target_compile_definitions(${name} PRIVATE "LOCAL_STAT=1") # switch off when do not collect local stat
    #    target_compile_definitions(${name} PRIVATE "SUP_STAT=1") # switch off when do not collect local stat
endfunction()

############### Parallel Latest Codes ################################
build_varying_parameters(pkt-varying "" "")
build_truss_algorithm(pkt-inter-shrink "PKT_INTER=1" pkt_merge.cpp)
build_truss_algorithm(demo "PKT_INTER=1" pkt_merge.cpp)
target_compile_definitions(demo PRIVATE "PAPER_FIGURE=1")
#target_compile_definitions(pkt-inter-shrink PRIVATE "SHRINK_EDGE_LIST=1")
#target_compile_definitions(pkt-inter-shrink PRIVATE "COMPACT_CSR=1")

set(PKT_MARK_SHRINK_FILES pkt_merge.cpp)

set(PKT_RADIX_HASH_SHRINK_FILES pkt_merge.cpp pkt_radix_map.cpp)
build_truss_algorithm(pkt-radix-shrink PKT_MARK=1 "${PKT_RADIX_HASH_SHRINK_FILES}")

################################################################

################## Experimental of Pre-Processing and Support Initialization ################################
build_eid_ours(pkt-eid-parallel "EID_TEST=1" pkt_merge.cpp)
build_eid_ours(pkt-eid-parallel-seq-scan "EID_TEST=1;SEQ_SCAN=1" pkt_merge.cpp)
build_eid_ours(pkt-eid-parallel-rev "EID_TEST=1;REVERSE=1" pkt_merge.cpp)
add_executable(pkt-varying-wpt exp_varying_parameters/main_eval_tc_varing_wpt.cpp)
target_link_libraries(pkt-varying-wpt PRIVATE graph-pre-processing)

build_eval_tc(pkt-eval-tc-baseline "" "")
build_eval_tc(pkt-eval-tc-dtc "DTC=1" "")
build_eval_tc(pkt-eval-tc-wp "WP=1" "")
################################################################

if (PLAYGROUND)
    add_subdirectory(playground)
endif ()
