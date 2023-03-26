#ifndef ADSORPTION_POTENTIALS_H
#define ADSORPTION_POTENTIALS_H

#include <cmath>
#include "helpers.h"
#include "data_classes.h"

double DRA(double z, double eps0, double z0, double beta);
double LEE(double z, double eps_k, double sigma_ff, Adsorbent adsorbent);
double STEELE(double z, double eps_k, double sigma_ff, Adsorbent adsorbent);
#endif