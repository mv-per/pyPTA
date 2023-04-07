cmake_minimum_required(VERSION 3.16)

set(This cPTA)

project(${This} VERSION 0.1.0 )

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(gtest_force_shared_crt ON)

enable_testing()

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)


add_subdirectory(googletest)
add_subdirectory(tests)

file (GLOB Sources "src/*.cpp")
file (GLOB OptimizationSources "src/optimization_algorithms/brent.cpp")
file (GLOB EoSSources "src/equations_of_state/*.cpp")

file (GLOB OptimizationHeaders "src/optimization_algorithms/brent.h")
file (GLOB EoSHeaders "src/equations_of_state/*.h")
file (GLOB Headers "src/*.h")


add_library(${This} SHARED ${Sources} ${OptimizationSources} ${EoSSources} ${Headers} ${OptimizationHeaders} ${EoSHeaders})
