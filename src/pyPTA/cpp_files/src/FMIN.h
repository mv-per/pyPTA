const double C_FMIN = 0.5*(3. - sqrt(5.0));


double brent_fmin(double (*fun)(double), double x, double tol)
{
	double a, b, d, e, xm, p, q, r, tol1, tol2, u, v, w, eps, X;
    double fu, fv, fw, fx;

	//Compute eps, the relative machine precision
	// eps = 1.1e-16;


	
	// printf("%E\n", tol);
    getboundaries(fun, x, &a, &b);

    //If manually determination of eps
    // eps = 1.;
	// define_eps:
	// 	eps = eps/2.0;
	// 	tol1 = 1.0 + eps;
	// 	if (tol1 > 1.0) goto define_eps;
    //     eps = sqrt(eps);

    // printf("%e\n", eps);

    // Initialization
    v = a + C_FMIN * (b-a);
    w = v;
    X = v;
    e = 0.;
    fx = fun(X);
    fv = fx;
    fw = fx;
    
	//Begin step
	begin_20:
        xm = 0.5 * (a + b);
        tol1 = eps*fabs(X) + tol/3.0;
        tol2 = 2.0*tol1;

    //Check stop criterion
    if (fabs(X-xm) <= (tol2 - 0.5*(b-a))) goto done;

    // Is golden search necessary?
    if (fabs(e) <= tol1) goto begin_40;

    // Fit parabola
    r = ( X- w) * (fx - fv);
    q = (X - v) * (fx - fw);
    p = (X - v) * q - (X - w) * r;
    q = 2.0 * (q - r);
    if (q > 0.0) p = -p;
    q = fabs(q);
    r = e;
    e = d;

    // Is parabola acceptable
	begin_30:
		if (fabs(p) >= fabs(0.5 * q * r)) goto begin_40;
        if (p <= q * (a - X)) goto begin_40;
        if (p >= q * (b - X)) goto begin_40;

    // A parabolic interpolation step
    d = p/q;
    u = X + d;

    // fun must not be evaluated too close to a or b
    if ((u-a) < tol2) d = sign(tol1, xm-X);
    if ((b-u) < tol2) d = sign(tol1, xm-X);
    goto begin_50;

	begin_40:
        if (X >= xm) e = a - X;
        if (X < xm) e = b - X;
        d = C_FMIN * e;

    begin_50:
        if (fabs(d) >= tol1) u = X + d;
        if (fabs(d) < tol1) u = X + sign(tol1, d);
        fu = fun(u);
    
    // update a, b, v, w, and X
    if (fu > fx) goto begin_60;
    if (u >= X) a = X;
    if (u < X) b = X;
    v = w;
    fv = fw;
    w = X;
    X = u;
    fx = fu;
    goto begin_20;

    begin_60:
        if (u < X) a = u;
        if (u >= X) b = u;
        if (fu <= fw) goto begin_70;
        if (w == X) goto begin_70;
        if (fu <= fv) goto begin_80;
        if (v == X) goto begin_80;
        if (v == w) goto begin_80;
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