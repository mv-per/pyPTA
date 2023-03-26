cmake_minimum_required(VERSION 3.25.0)

set(This pyPTATests)

set(Sources
    test_pta_mixture.cpp
    test_pta_pure.cpp
    test_logic.cpp
    test_loadings.cpp
)

add_executable(${This} ${Sources})
target_link_libraries(${This} PUBLIC
    gtest_main
    pyPTA
)

add_test(
    NAME ${This}
    COMMAND ${This}
)