#include "helpers.h"

/**
 *
 * Return the proper deviation function from the deviation type
 *
 * @param deviation_type Deviation type.
 * @return Deviation function.
 */
deviation_caller get_deviation_function(std::string deviation_type)
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

std::vector<double> linespace(double start, double ed, std::size_t num)
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
