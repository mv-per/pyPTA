

#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>
#include "equations_of_state/eos.h"
#include "optimization_algorithms/brent.h"
#include "data_classes.h"

class PTASolver
{
private:
    double BulkFugacity, PotentialEnergy, Temperature;
    call_mono_eos EquationOfState;
    call_mix_eos MixEquationOfState;
    double OptimizedPressure;

    /// Default tolerance
    double DEFAULT_TOL = 1e-12;

    double equilibrium(double p_,
                       double f_eps);

public:
    PTASolver(double bulk_fugacity, double temperature, call_mono_eos eos);
    PTASolver(double bulk_fugacity, double temperature, call_mix_eos eos);
    double findOptimizedPressure(double InitialEstimate, double f_eps);
};
#endif