

#include "../global_helper.h"
#include <iostream>
#include <cassert>

int main(void)
{

    Deviation absolute_deviation = get_deviation_function("absolute");
    Deviation relative_deviation = get_deviation_function("relative");
    Deviation absolute_relative_deviation = get_deviation_function("absolute_relative");

    assert(absolute_deviation(0.1, 1) == 0.9);
    std::cout << "Expression valid....\n"
              << std::endl;

    assert(relative_deviation(0.1, 1) == -9.);
    std::cout << "Expression valid....\n"
              << std::endl;

    assert(absolute_relative_deviation(0.1, 1) == 9.);
    std::cout << "Expression valid....\n"
              << std::endl;

    return 0;
}