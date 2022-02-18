double BISECTION(double (*fun)(double), double x, double tol)
{
    double xn = x;
    double a, b;
    double fa;
    double c, fc;
 
	// find a and b
    getboundaries(fun, x, &a, &b);

    // printf("a: %f, b: %f\n", a, b);
	//Initialization
    size_t k = 0;
    do {
        c = 0.5*(a+b);
        fc = fun(c);
        fa = fun(a);

        if ((fc == 0) || ((b-a)/2.0 <= tol )) break;

        if ((fc > 0.) && (fa >0)) a = c;
        else b = c;

        k++;
        // printf("%d\n", k);
    } while (k <= 150);

    return c;
}