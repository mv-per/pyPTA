#include "pta_solver.h"

#include <chrono>
#include <condition_variable>
#include <ctime>
#include <functional> // for _1, _2
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>

PTASolver::PTASolver(double bulk_fugacity, double temperature, call_mono_eos eos)
{
    this->EquationOfState = eos;
    this->BulkFugacity = bulk_fugacity;
    this->Temperature = temperature;
}

PTASolver::PTASolver(double bulk_fugacity, double temperature, call_mix_eos eos)
{
    this->MixEquationOfState = eos;
    this->BulkFugacity = bulk_fugacity;
    this->Temperature = temperature;
}

double PTASolver::equilibrium(double P,
                              double f_eps)
{
    if (std::isinf(P) == 1 || std::isnan(P) == 1)
    {
        return 10000.;
    }
    struct mono_eos adsorbed = this->EquationOfState(P, this->Temperature);
    return this->BulkFugacity * f_eps - adsorbed.fug;
}

double PTASolver::findOptimizedPressure(double InitialPressureEstimate, double f_eps)
{

    auto equilibrium = [&](double x)
    {
        return this->equilibrium(x, f_eps) * 1000;
    };

    return brent_zeroin(equilibrium, InitialPressureEstimate, this->DEFAULT_TOL);
}