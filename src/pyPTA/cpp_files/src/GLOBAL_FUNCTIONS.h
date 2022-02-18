#pragma once
#include <iostream>
using namespace std;
#include <vector>


/* vetor linear */ vector<double> linespace(double start, double ed, size_t num);
/* size_tegração */ double SIMPS(vector<double> vector, double step, size_t maxlength);
/* Minimização */ double brent_zeroin(double (*fun)(double), double x, double tol);
/* função sign */ double sign(double X1, double X2);

double minVal(vector<double> valVec)
{
	double MIN = valVec[0];
	for (size_t i = 0; i < valVec.size(); i++)
	{
		if (valVec[i] < MIN)
		{
			MIN = valVec[i];
		}
	}
	return MIN;
}

vector<double> linespace(double start, double ed, size_t num)
{
	// catch rarely, throw often
	if (num < 2)
	{
		throw new exception();
	}
	size_t partitions = num - 1;
	vector<double> pts;
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

double SIMPS(vector<double> vector, double step, size_t maxlength)
{
	size_t i;
	double so = 0;
	double se = 0;

	for (i = 0; i < maxlength; i++)
	{
		if (i % 2 == 1)
		{
			if (vector[i] != vector[i])
			{
				vector[i] = vector[i - 1];
			}
			so += vector[i];
		}
		else
		{
			if (vector[i] != vector[i])
			{
				vector[i] = vector[i - 1];
			}
			se += vector[i];
		}
	}
	double result = step / 3.0 * (vector[0] + vector[maxlength - 1] + 4.0 * so + 2.0 * se);
	return result;
}



double sign(double X1, double X2)
{
	if (X2 >= 0)
		return fabs(X1);
	else
		return -fabs(X1);
}

double brent_zeroin(double (*fun)(double), double x, double tol)
{
	double a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s, eps, fx, dx;
	//Compute eps, the relative machine precision
	eps = 1.1e-16;

	/* If manually determination of eps
	define_eps:
		eps = eps/2.0;
		tol1 = 1.0 + eps;
		if (tol1 > 1.0) goto define_eps;
	std::cout << eps << std::endl;
	*/
	fx = fun(x);
	// find a and b
	/*
	if (fx == 0){
		fb = fx;
		b = x;
	    goto done;
	}
	*/

	if (x != 0.0)
		dx = x / 50.0;
	else
		dx = 1.0 / 50.0;

	double sqrt_2 = 1.414213562373095048801;
	//assign initial values
	a = x;
	fa = fx;
	b = x;
	fb = fx;

	// find change of sign - copiado do matlab
	while ((fa > 0) == (fb > 0))
	{
		dx = sqrt_2 * dx;
		a = x - dx;
		fa = fun(a);

		if ((fa > 0) != (fb > 0))
			break;

		b = x + dx;
		fb = fun(b);
	}

	//Initialization

	fa = fun(a);
	fb = fun(b);

//Begin step
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
			cout << array << endl;
		}
	}
	return array;
}