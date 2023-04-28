#include "bisection.h"

void getboundaries(std::function<double(double)> fun, double &x, double *a, double *b)
{
    double fa, fb, fx, dx;
    // assign initial values
    fx = fun(x);
    *a = x;
    fa = fx;
    *b = x;
    fb = fx;

    // find a and b
    if (x != 0.)
        dx = x / 50.0;
    else
        dx = 1.0 / 50.0;

    std::size_t k = 0;
// find change of sign - copiado do matlab
get_a_and_b:
    while ((fa > 0.) == (fb > 0.))
    {
        dx = 1.414213562373095048801 * dx;
        *a = x - dx;
        fa = fun(*a);
        if ((fa > 0.) != (fb > 0.))
            break;
        *b = x + dx;
        fb = fun(*b);
    }

    if ((std::fabs(fa) > .5 || fabs(fb) > .5) && k < 10)
    {
        dx /= 2.;
        k++;
        goto get_a_and_b;
    }
}

double bisection(std::function<double(double)> fun, double x, double tol)
{
    double xn = x;
    double a, b;
    double fa;
    double c, fc;

    // find a and b
    getboundaries(fun, x, &a, &b);

    // printf("a: %f, b: %f\n", a, b);
    // Initialization
    std::size_t k = 0;
    do
    {
        c = 0.5 * (a + b);
        fc = fun(c);
        fa = fun(a);

        if ((fc == 0) || ((b - a) / 2.0 <= tol))
            break;

        if ((fc > 0.) && (fa > 0))
            a = c;
        else
            b = c;

        k++;
        // printf("%d\n", k);
    } while (k <= 150);

    return c;
}
