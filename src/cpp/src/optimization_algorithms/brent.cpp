
#include "brent.h"
#include <iostream>

/**
 * Returns the first value in the same unit as de second one
 *
 * @param value1 Value to change sign.
 * @param value2 Value to find sign.
 * @return value1 in the same unit as value2
 */
double sign(double value1, double value2)
{
    if (value2 >= 0)
        return fabs(value1);
    else
        return -fabs(value1);
}

/**
 * Perform the Richard Brent's improvements to Dekker's zeroin algorithm
 *
 * @param fun Pointer to Function that receives a double and returns a double (function to be minimized).
 * @param x Initial estimate.
 * @param tol Minimum `TOL` .
 * @return Value of the value with minimal :fun: value.
 */
double brent_zeroin(std::function<double(double)> fun, double x, double tol)
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
        dx = 1.414213562373095048801 * dx;
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
    return b;
}
