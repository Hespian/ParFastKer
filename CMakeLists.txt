cmake_minimum_required(VERSION 2.8.7)
project(fast-reductions)

# Inlcude for some machines
include_directories(/usr/include/mpi)
 
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Include directories for executable
include_directories(extern/KaHIP)
include_directories(extern/KaHIP/lib)
include_directories(extern/KaHIP/lib/tools)
include_directories(extern/KaHIP/lib/io)

include_directories(extern/KaHIP/interface)

#enable c++11
add_definitions(-std=c++11)

set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} -fopenmp")
add_definitions(${OpenMP_CXX_FLAGS})
add_definitions(${CMAKE_CXX_FLAGS_RELEASE})

# Add the test executable
set(TEST_SOURCES ${SOURCES})
list(REMOVE_ITEM TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/benchmark.cpp")
add_executable(test_correctness ${TEST_SOURCES})

# Add the benchmark executable
set(BENCHMARK_SOURCES ${SOURCES})
list(REMOVE_ITEM BENCHMARK_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/test_correctness.cpp")
add_executable(benchmark ${BENCHMARK_SOURCES})

# Additional include directories for libkahip
include_directories(extern/KaHIP/lib/partition)
include_directories(extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)

# Build the KaHIP library
file(GLOB_RECURSE lib_kahip_files "extern/KaHIP/*.cpp")
add_library(libkahip STATIC ${lib_kahip_files})

# Link libraries
target_link_libraries(test_correctness argtable2)
target_link_libraries(test_correctness libkahip)
target_link_libraries(test_correctness gomp)
target_link_libraries(benchmark argtable2)
target_link_libraries(benchmark libkahip)
target_link_libraries(benchmark gomp)