
#ifndef GLOBAL_HELPER_H
#define GLOBAL_HELPER_H

#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <stdexcept>
#include "constants.h"

/* Deviation caller type */
typedef std::function<double(double, double)> Deviation;

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
 * Perform the Richard Brent's improvements to Dekker's zeroin algorithm
 *
 * @param fun Pointer to Function that receives a double and returns a double (function to be minimized).
 * @param x Initial estimate.
 * @param tol Minimum `TOL` .
 * @return Value of the value with minimal :fun: value.
 */
double brent_zeroin(double (*fun)(double), double x, double tol);

/**
 * Returns the first value in the same unit as de second one
 *
 * @param value1 Value to change sign.
 * @param value2 Value to find sign.
 * @return value1 in the same unit as value2
 */
double sign(double value1, double value2);

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

//-------------------------

/**
 *
 * Return the proper deviation function from the deviation type
 *
 * @param deviation_type Deviation type.
 * @return Deviation function.
 */
std::function<double(double, double)> get_deviation_function(std::string deviation_type)
{
	if (deviation_type == "absolute")
	{
		return [](double A, double B)
		{ return absolute_deviation(A, B); };
	}
	else if (deviation_type == "relative")
	{
		return [](double A, double B)
		{ return relative_deviation(A, B); };
	}
	else if (deviation_type == "absolute_relative")
	{
		return [](double A, double B)
		{ return absolute_relative_deviation(A, B); };
	}
	else
	{
		throw std::invalid_argument("Invalid Deviation Type");
	}
}

double
absolute_relative_deviation(double A, double B)
{
	return std::fabs((A - B) / A);
}
double
relative_deviation(double A, double B)
{
	return (A - B) / A;
}

double
absolute_deviation(double A, double B)
{
	return std::fabs(A - B);
}

double min_vec_val(std::vector<double> values)
{
	double MIN = values[0];
	for (size_t i = 0; i < values.size(); i++)
	{
		if (values[i] < MIN)
		{
			MIN = values[i];
		}
	}
	return MIN;
}

std::vector<double> linespace(double start, double ed, size_t num)
{
	// catch rarely, throw often
	if (num < 2)
	{

		throw std::invalid_argument("Number of points must be bigger or equal than 2");
	}
	size_t partitions = num - 1;
	std::vector<double> pts;
	// length of each segment
	double length = (ed - start) / double(partitions);
	// first, not to change
	pts.push_back(start);
	for (size_t i = 1; i < num - 1; i++)
	{
		pts.push_back(start + i * length);
	}
	// last, not to change
	pts.push_back(ed);
	return pts;
}

double simpson_rule(std::vector<double> integral_vector, double step)
{
	std::size_t maxlength = integral_vector.size();
	std::size_t i;
	double so = 0;
	double se = 0;

	for (i = 0; i < maxlength; i++)
	{
		if (i % 2 == 1)
		{
			if (integral_vector[i] != integral_vector[i])
			{
				integral_vector[i] = integral_vector[i - 1];
			}
			so += integral_vector[i];
		}
		else
		{
			if (integral_vector[i] != integral_vector[i])
			{
				integral_vector[i] = integral_vector[i - 1];
			}
			se += integral_vector[i];
		}
	}
	double result = step / 3.0 * (integral_vector[0] + integral_vector[maxlength - 1] + 4.0 * so + 2.0 * se);
	return result;
}

double sign(double value1, double value2)
{
	if (value2 >= 0)
		return fabs(value1);
	else
		return -fabs(value1);
}

double brent_zeroin(double (*fun)(double), double x, double tol)
{
	double a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s, eps, fx, dx;

	eps = 1.1e-16;

	fx = fun(x);

	if (x != 0.0)
		dx = x / 50.0;
	else
		dx = 1.0 / 50.0;

	// assign initial values
	a = x;
	fa = fx;
	b = x;
	fb = fx;

	// Finds a sign change
	while ((fa > 0) == (fb > 0))
	{
		dx = SQRT_2 * dx;
		a = x - dx;
		fa = fun(a);

		if ((fa > 0) != (fb > 0))
			break;

		b = x + dx;
		fb = fun(b);
	}

	// Initialization

	fa = fun(a);
	fb = fun(b);

// Begin step
begin_20:
	c = a;
	fc = fa;
	d = b - a;
	e = d;
begin_30:
	if (fabs(fc) >= fabs(fb))
		goto convergence_test;
	a = b;
	b = c;
	c = a;
	fa = fb;
	fb = fc;
	fc = fa;

convergence_test:
	tol1 = 2.0 * eps * fabs(b) + 0.5 * tol;
	xm = 0.5 * (c - b);
	if (fabs(xm) <= tol1 || fb == 0.0)
		goto done;

	// Is bisection necessary
	if (fabs(e) < tol1)
		goto bisection;
	if (fabs(fa) <= fabs(fb))
		goto bisection;

	// Is quadratic size_terpolation possible
	if (a != c)
		goto inverse_quadratic_size_terpolation;

	s = fb / fa;
	p = 2.0 * xm * s;
	q = 1.0 - s;

	goto adjust_signs;

inverse_quadratic_size_terpolation:
	q = fa / fc;
	r = fb / fc;
	s = fb / fa;
	p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
	q = (q - 1.0) * (r - 1.0) * (s - 1.0);

adjust_signs:
	if (p > 0.0)
	{
		q = -q;
	}
	p = fabs(p);

	// Is size_terpolation acceptable
	if (2.0 * p >= (3.0 * xm * q - fabs(tol1 * q)) || p >= fabs(0.5 * e * q))
		goto bisection;
	e = d;
	d = p / q;
	goto complete_step;

bisection:
	d = xm;
	e = d;

complete_step:
	a = b;
	fa = fb;
	if (fabs(d) > tol1)
		b = b + d;
	else
		b = b + sign(tol1, xm);
	fb = fun(b);
	if ((fb * (fc / fabs(fc))) > 0.0)
		goto begin_20;
	goto begin_30;

done:
	double zeroin = b;
	return zeroin;
}

template <class T>
T removenegative(T oldArray)
{
	T array;
	for (size_t i = 0; i < oldArray.size(); i++)
	{
		if (oldArray[i] > 0.0)
		{
			array[i] = oldArray[i];
		}
	}
	return array;
}

#endif