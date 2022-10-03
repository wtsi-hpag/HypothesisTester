#pragma once

#include <vector>
#include <math.h>
std::vector<double> means;
std::vector<double> errors;

double func(const std::vector<double> & x)
{
	double v = 1;
	double invSqrt = 1.0/sqrt(2*M_PI);
	for (int i = 0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];
		v*=invSqrt/errors[i] * exp(-0.5*d*d);
	}
	return v;
}
double logfunc(const std::vector<double> & x)
{
	double v = 0;
	double log_invSqrt = -0.5*log(2*M_PI); //1.0/sqrt(2*M_PI);
	for (int i = 0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];
		v+= log_invSqrt - log(errors[i]) -0.5*d*d;
	}
	return v;
}
double ddlogfunc(const std::vector<double> & x, int i, int j)
{
	if (i!=j)
	{
		return 0;
	}
	else
	{
		return -1.0/(errors[i] * errors[i]);
	}
}

void dlogfunc(std::vector<double> & grad, const std::vector<double> & x)
{
	for (int i = 0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];

		grad[i] = -d/(errors[i]);
	}
}
