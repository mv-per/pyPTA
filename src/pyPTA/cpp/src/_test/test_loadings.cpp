

#include <gtest/gtest.h>
#include "../helpers.h"

TEST(test_loadings, DemostrateDeviations)
{
    deviation_caller absolute_deviation = get_deviation_function("absolute");
    deviation_caller relative_deviation = get_deviation_function("relative");
    deviation_caller absolute_relative_deviation = get_deviation_function("absolute_relative");

    ASSERT_TRUE(absolute_deviation(0.1, 1) == 0.9);
    ASSERT_TRUE(relative_deviation(0.1, 1) == -9.);
    ASSERT_TRUE(absolute_relative_deviation(0.1, 1) == 9.);
}
