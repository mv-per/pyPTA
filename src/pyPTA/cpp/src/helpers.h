#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <stdexcept>

/* Deviation caller type */
typedef std::function<double(double, double)> deviation_caller;

deviation_caller get_deviation_function(std::string deviation_type);

/**
 * Creates a linear spaced vector from start to end with `num` steps
 *
 * @param start Initial value.
 * @param end Final Value
 * @param num Number of steps
 * @return Linear Spaced Vector
 */
std::vector<double> linespace(double start, double end, std::size_t num);

/**
 * Return the quadratic interpolation by the Simpson's 1/3 rule
 *
 * @param integral_vector Vector containing integral values.
 * @param step Step size of the vector.
 * @return Integrated value.
 */
double simpson_rule(std::vector<double> integral_vector, double step);

/**
 *
 * Finds the smallest value of a vector
 *
 * @param values Vector of double values.
 * @return Smallest value.
 */
double min_vec_val(std::vector<double> valVec);

/**
 *
 * Finds the absolute deviation from A and B
 *
 * @param A Value of A.
 * @param B Value of B.
 * @return Absolute deviation.
 */
double absolute_deviation(double A, double B);

/**
 *
 * Finds the relative deviation from A and B
 *
 * @param A Value of A.
 * @param B Value of B.
 * @return Absolute deviation.
 */

double relative_deviation(double A, double B);

/**
 *
 * Finds the absolute relative deviation from A and B
 *
 * @param A Value of A.
 * @param B Value of B.
 * @return Absolute deviation.
 */
double absolute_relative_deviation(double A, double B);

#endif