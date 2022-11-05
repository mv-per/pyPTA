#include <iostream>
#include <functional>
#include "../pr77.h"

typedef std::function<mono_eos(double, double)> call_eos;

std::function<mono_eos(double, double)> get_equation_of_state(std::vector<double> eos_params)
{
    return [=](double P, double T)
    {
        return pr77().get_mono_fluid_properties(P, T, eos_params[0], eos_params[1], eos_params[2]);
    };
}

int main()
{
    call_eos eos_pointer = get_equation_of_state({73e5, 303, 0.222});

    cout << eos_pointer(1000000, 303).fug << endl;

    return 0;
}
