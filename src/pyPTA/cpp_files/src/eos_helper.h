

double R = 8.314462618; // m3 Pa /K / mol

struct Fluid
{
	double Pc;
	double Tc;
	double w;
};

struct mono_eos
{
	double fug, dens, phi, Z;
};

struct mix_eos
{
	vector<double> fug;
	double dens;
	vector<double> phi;
	double Z;
};

vector<double> find_z(double A, double B);
double gx(double X, double a0, double a1, double a2);
double dgx(double X, double a1, double a2);
double d2gx(double X, double a2);
double minvalue(double num1, double num2, double num3);
double maxvalue(double num1, double num2, double num3);

vector<double> find_z(double A, double B)
{

	double a2, a1, a0, amax, rr, xinf1, xold, xnew, x1, x2, x3, gx0, gx1, gx2, a, b, c;
	std::size_t iter;
	a2 = B - 1.0;
	a1 = A - 3.0 * B * B - 2.0 * B;
	a0 = -A * B + B * B + B * B * B;
	amax = maxvalue(fabs(a2), fabs(a1), fabs(a0));
	rr = 1.0 + amax;
	xinf1 = -a2 / 3.0;

	if (gx(xinf1, a0, a1, a2) > 0.0)
	{
		xold = -rr;
	}
	else
	{
		xold = rr;
	}

	iter = 0;
	while (iter < 50)
	{
		gx0 = gx(xold, a0, a1, a2);
		gx1 = dgx(xold, a1, a2);
		gx2 = d2gx(xold, a2);
		xnew = xold - gx0 * gx1 / (gx1 * gx1 - 1.0 / 2.0 * gx0 * gx2); // kepler method

		if ((fabs(xnew - xold) < 1.0e-10) || (gx0 == 0))
		{
			break;
		}
		xold = xnew;
		iter++;
	}

	x1 = xnew;
	a = 1.0;
	b = a * x1 + a2;
	c = b * x1 + a1;

	double discrim = b * b - 4.0 * a * c;
	if (discrim < 0.0)
	{
		x2 = x1;
		x3 = x1;
	}
	else
	{
		x2 = (-b + sqrt(discrim)) / (2.0 * a);
		x3 = (-b - sqrt(discrim)) / (2.0 * a);
	}

	vector<double> results(3, 0.0);
	results = {x1, x2, x3};

	return results;
}

double gx(double X, double a0, double a1, double a2)
{
	return X * X * X + a2 * X * X + a1 * X + a0;
}

double dgx(double X, double a1, double a2)
{
	return 3.0 * X * X + 2.0 * a2 * X + a1;
}

double d2gx(double X, double a2)
{
	return 6.0 * X + 2.0 * a2;
}

double minvalue(double num1, double num2, double num3)
{
	double result;
	if ((num1 < num2) && (num1 < num3) && (fabs(num1) == num1) && (fabs(num2) == num2) && (fabs(num3) == num3))
	{
		result = num1;
	}
	else if ((num2 < num1) && (num2 < num3) && (fabs(num1) == num1) && (fabs(num2) == num2) && (fabs(num3) == num3))
	{
		result = num2;
	}
	else
	{
		result = num3;
	}
	return result;
}

double maxvalue(double num1, double num2, double num3)
{
	double result;
	if ((num1 > num2) && (num1 > num3))
	{
		result = num1;
	}
	else if ((num2 > num1) && (num2 > num3))
	{
		result = num2;
	}
	else
	{
		result = num3;
	}
	return result;
}