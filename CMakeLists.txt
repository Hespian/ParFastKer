cmake_minimum_required(VERSION 3.5.2)
project(fast-reductions)
 
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Include directories for executable
include_directories(extern/KaHIP)
include_directories(extern/KaHIP/lib)
include_directories(extern/KaHIP/lib/tools)
include_directories(extern/KaHIP/lib/io)

include_directories(extern/KaHIP/interface)

#enable c++11
add_definitions(-std=c++11)
 
add_definitions(${OpenMP_CXX_FLAGS})
add_definitions(${CMAKE_CXX_FLAGS_RELEASE})

# Add the executable
add_executable(fast-reductions ${SOURCES})

# Additional include directories for libkahip
include_directories(extern/KaHIP/lib/partition)
include_directories(extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)

# Build the KaHIP library
file(GLOB_RECURSE lib_kahip_files "extern/KaHIP/*.cpp")
add_library(libkahip STATIC ${lib_kahip_files})

# Link libraries
target_link_libraries(fast-reductions argtable2)
target_link_libraries(fast-reductions libkahip)
target_link_libraries(fast-reductions gomp)