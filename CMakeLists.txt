cmake_minimum_required(VERSION 3.4...3.18)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


#  Project
set(This pyPTA)

project(${This})

include_directories("${CMAKE_SOURCE_DIR}/src/cpp/src")
include_directories("${CMAKE_SOURCE_DIR}/src/cpp/src/optimization_algorithms")
include_directories("${CMAKE_SOURCE_DIR}/src/cpp/src/equations_of_state")

file (GLOB Sources "src/cpp/src/*.cpp")
file (GLOB Headers "src/cpp/src/*.h")
file (GLOB OptimizationSources "src/cpp/src/optimization_algorithms/*.cpp")
file (GLOB OptimizationHeaders "src/cpp/src/*.h")
file (GLOB EOSSources "src/cpp/src/equations_of_state/*.cpp")
file (GLOB EOSHeaders "src/cpp/src/equations_of_state/*.h")

find_package(Python COMPONENTS Interpreter Development)
add_subdirectory(/pybind11)
# find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(
    ${This} 
    ${Sources} 
    ${OptimizationSources}
    ${EOSSources}
    ${Headers} 
    ${OptimizationHeaders} 
    ${EOSHeaders} 
    src/module.cpp)

target_link_libraries(${This} PUBLIC)


# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(${This}
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})