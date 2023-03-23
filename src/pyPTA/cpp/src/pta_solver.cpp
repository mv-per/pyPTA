#include "pta_solver.h"

PTASolver::PTASolver(double bulk_fugacity, double temperature, call_mono_eos eos)
{
    this->EquationOfState = eos;
    this->BulkFugacity = bulk_fugacity;
    this->Temperature = temperature;
}

double PTASolver::find_pz(double InitialEstimate, double f_eps)
{
    double OptimizedPressure = 0.0;

    static auto get_equilibrium = [&](double x)
    {
        return this->equilibrium(x, this->BulkFugacity, f_eps, this->Temperature, this->EquationOfState);
    };
    double (*equilibrium)(double) = [](double p)
    { return get_equilibrium(p); };
    OptimizedPressure = brent_zeroin(equilibrium, InitialEstimate, this->DEFAULT_TOL);

    return OptimizedPressure;
}

double PTASolver::equilibrium(double p_,
                              double bulk_fugacity,
                              double f_eps,
                              double T,
                              call_mono_eos eos)
{
    if (std::isinf(p_) == 1 || std::isnan(p_) == 1)
    {
        return 10000.;
    }
    struct mono_eos adsorbed = eos(p_, T);
    return bulk_fugacity * f_eps - adsorbed.fug;
}
