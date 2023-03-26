cmake_minimum_required(VERSION 3.26.0)

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

# Headers
set(PROJECT_SOURCE_DIR src)



add_subdirectory(googletest)
add_subdirectory(${PROJECT_SOURCE_DIR}/_test)


set(Headers
    ${PROJECT_SOURCE_DIR}/helpers.h
    ${PROJECT_SOURCE_DIR}/data_classes.h
    ${PROJECT_SOURCE_DIR}/adsorption_potentials.h
    ${PROJECT_SOURCE_DIR}/pta_solver.h
    ${PROJECT_SOURCE_DIR}/pta_pure.h
    ${PROJECT_SOURCE_DIR}/pta_mixture.h
    ${PROJECT_SOURCE_DIR}/pta_helper.h
    ${PROJECT_SOURCE_DIR}/equations_of_state/eos.h
    ${PROJECT_SOURCE_DIR}/equations_of_state/eos_helper.h
    ${PROJECT_SOURCE_DIR}/equations_of_state/pr77.h
    ${PROJECT_SOURCE_DIR}/optimization_algorithms/brent.h
    ${PROJECT_SOURCE_DIR}/optimization_algorithms/nelder_mead.h
)

set(Sources
    ${PROJECT_SOURCE_DIR}/helpers.cpp
    ${PROJECT_SOURCE_DIR}/adsorption_potentials.cpp
    ${PROJECT_SOURCE_DIR}/pta_solver.cpp
    ${PROJECT_SOURCE_DIR}/pta_pure.cpp
    ${PROJECT_SOURCE_DIR}/pta_mixture.cpp
    ${PROJECT_SOURCE_DIR}/pta_helper.cpp
    ${PROJECT_SOURCE_DIR}/equations_of_state/eos_helper.cpp
    ${PROJECT_SOURCE_DIR}/equations_of_state/pr77.cpp
    ${PROJECT_SOURCE_DIR}/optimization_algorithms/brent.cpp
)

# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${Sources})

add_library(${This} STATIC ${Sources} ${Headers})
# add_library(${This} SHARED ${Sources})


# # Include directories
# target_include_directories(automobile PRIVATE include/automobile_bits)

# Install
# install(TARGETS ${This} DESTINATION lib)

# Install the headers
# install(FILES ${Headers} DESTINATION include)

# Create base directory
# install(DIRECTORY ${Headers} DESTINATION include)