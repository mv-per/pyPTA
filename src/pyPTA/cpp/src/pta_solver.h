

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
    double BulkFugacity, PotentialEnergy, Temperature;
    call_mono_eos EquationOfState;
    /// Default tolerance
    double DEFAULT_TOL = 1e-12;

public:
    PTASolver(double bulk_fugacity, double temperature, call_mono_eos eos);
    double find_pz(double InitialEstimate, double f_eps);
    double equilibrium(double p_,
                       double bulk_fugacity,
                       double f_eps,
                       double T,
                       call_mono_eos eos);
};
#endif