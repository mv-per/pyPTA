#ifndef BISECTION_H
#define BISECTION_H

#include <functional>
#include <cmath>

double bisection(std::function<double(double)> fun, double x, double tol);

#endif