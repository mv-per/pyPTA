#include "fmin.h"

const double C_FMIN = 0.5 * (3. - sqrt(5.0));

std::tuple<double, double> fmin_getboundaries(std::function<double(double)> fun, double &x)
{
    double fa, fb, fx, dx, a, b;
    // assign initial values
    fx = fun(x);
    a = x;
    fa = fx;
    b = x;
    fb = fx;

    // find a and b
    if (x != 0.)
        dx = x / 50.0;
    else
        dx = 1.0 / 50.0;

    size_t k = 0;
// find change of sign - copiado do matlab
get_a_and_b:
    while ((fa > 0.) == (fb > 0.))
    {
        dx = 1.414213562373095048801 * dx;
        a = x - dx;
        fa = fun(a);
        if ((fa > 0.) != (fb > 0.))
            break;
        b = x + dx;
        fb = fun(b);
    }

    if ((fabs(fa) > .5 || fabs(fb) > .5) && k < 10)
    {
        dx /= 2.;
        k++;
        goto get_a_and_b;
    }

    return std::make_tuple(a, b);
}

/**
 * Returns the first value in the same unit as de second one
 *
 * @param value1 Value to change sign.
 * @param value2 Value to find sign.
 * @return value1 in the same unit as value2
 */
double fmin_sign(double value1, double value2)
{
    if (value2 >= 0)
        return fabs(value1);
    else
        return -fabs(value1);
}

double brent_fmin(std::function<double(double)> fun, double x, double tol)
{
    double a, b, d, e, xm, p, q, r, tol1, tol2, u, v, w, eps, X, dx;
    double fu, fv, fw, fx;

    // Compute eps, the relative machine precision
    //  eps = 1.1e-16;

    std::tie(a, b) = fmin_getboundaries(fun, x);

    // If manually determination of eps
    //  eps = 1.;
    //  define_eps:
    //  	eps = eps/2.0;
    //  	tol1 = 1.0 + eps;
    //  	if (tol1 > 1.0) goto define_eps;
    //      eps = sqrt(eps);

    // printf("%e\n", eps);

    // Initialization
    v = a + C_FMIN * (b - a);
    w = v;
    X = v;
    e = 0.;
    fx = fun(X);
    fv = fx;
    fw = fx;

// Begin step
begin_20:
    xm = 0.5 * (a + b);
    tol1 = eps * fabs(X) + tol / 3.0;
    tol2 = 2.0 * tol1;

    // Check stop criterion
    if (fabs(X - xm) <= (tol2 - 0.5 * (b - a)))
        goto done;

    // Is golden search necessary?
    if (fabs(e) <= tol1)
        goto begin_40;

    // Fit parabola
    r = (X - w) * (fx - fv);
    q = (X - v) * (fx - fw);
    p = (X - v) * q - (X - w) * r;
    q = 2.0 * (q - r);
    if (q > 0.0)
        p = -p;
    q = fabs(q);
    r = e;
    e = d;

// Is parabola acceptable
begin_30:
    if (fabs(p) >= fabs(0.5 * q * r))
        goto begin_40;
    if (p <= q * (a - X))
        goto begin_40;
    if (p >= q * (b - X))
        goto begin_40;

    // A parabolic interpolation step
    d = p / q;
    u = X + d;

    // fun must not be evaluated too close to a or b
    if ((u - a) < tol2)
        d = fmin_sign(tol1, xm - X);
    if ((b - u) < tol2)
        d = fmin_sign(tol1, xm - X);
    goto begin_50;

begin_40:
    if (X >= xm)
        e = a - X;
    if (X < xm)
        e = b - X;
    d = C_FMIN * e;

begin_50:
    if (fabs(d) >= tol1)
        u = X + d;
    if (fabs(d) < tol1)
        u = X + fmin_sign(tol1, d);
    fu = fun(u);

    // update a, b, v, w, and X
    if (fu > fx)
        goto begin_60;
    if (u >= X)
        a = X;
    if (u < X)
        b = X;
    v = w;
    fv = fw;
    w = X;
    X = u;
    fx = fu;
    goto begin_20;

begin_60:
    if (u < X)
        a = u;
    if (u >= X)
        b = u;
    if (fu <= fw)
        goto begin_70;
    if (w == X)
        goto begin_70;
    if (fu <= fv)
        goto begin_80;
    if (v == X)
        goto begin_80;
    if (v == w)
        goto begin_80;
    goto begin_20;

begin_70:
    v = w;
    fv = fw;
    w = u;
    fw = fu;
    goto begin_20;

begin_80:
    v = u;
    fv = fu;
    goto begin_20;

done:
    // printf("f(x) = %.8f \n", fabs(xm));
    double fmin = b;
    return fmin;
}