cmake_minimum_required(VERSION 2.8)
project(TC)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(SHARED FALSE CACHE BOOL "build a shared library")

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-D__OPENMP__ -DUSE_LOG)
add_compile_options(-g)

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

if (MSVC)
    set(TC_INSTALL FALSE)
else ()
    set(TC_INSTALL TRUE)
endif ()

# Add GK's standard cmake settings
include(./conf/gkbuild.cmake)

# Source files
file(GLOB tc_sources *.c)

# Include directories
include_directories(.)
#include_directories(~/local/include)

include_directories(../2017-hepc-K-Truss/GKlib)
add_subdirectory(../2017-hepc-K-Truss/GKlib gklib)
if (NOT GKLIB_PATH)
    message("NOT")
    set(GKLIB_PATH ../2017-hepc-K-Truss/GKlib)
endif ()

if (GKLIB_PATH)
    include_directories(${GKLIB_PATH}/include)
endif (GKLIB_PATH)

# Link directories 
#link_directories(~/local/lib)
#if(GKLIB_PATH)
#  link_directories(${GKLIB_PATH}/lib)
#endif(GKLIB_PATH)

# Build
set(PKGNAME tc-karypis)
set(UTIL_FILES yche/log.cpp yche/yche_helpers.cpp)
add_executable(${PKGNAME} ${tc_sources} ${UTIL_FILES})

foreach (prog ${PKGNAME})
    target_link_libraries(${prog} GKlib m)
endforeach (prog)

# Install
if (TC_INSTALL)
    install(TARGETS ${PKGNAME}
            RUNTIME DESTINATION bin)
endif ()
