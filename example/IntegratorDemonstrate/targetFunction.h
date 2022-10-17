#pragma once

#include <vector>
#include <math.h>
std::vector<double> means;
std::vector<double> errors;

double peak = 1;

inline double gaussFunc(const std::vector<double> & x)
{
	double v = 1;
	double invSqrt = 1.0/sqrt(2*M_PI);
	for (int i = 0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];
		v*=invSqrt/errors[i] * exp(-0.5*d*d);
	}
	return peak*v;
} 




double func(const std::vector<double> & x)
{
	// double v = 1;
	// double invSqrt = 1.0/sqrt(2*M_PI);
	// for (int i = 0; i < x.size(); ++i)
	// {
	// 	double d = (x[i] - means[i])/errors[i];
	// 	v*=invSqrt/errors[i] * exp(-0.5*d*d);
	// }
	// return peak*v;
	double v = 1;
	double inv = 1.0/M_PI;
	for (int i = 0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];
		v*= inv/errors[i] * 1.0/(1 + d*d);
	}
	return v;
}
double logfunc(const std::vector<double> & x)
{
	// double v = 0;
	// double log_invSqrt = -0.5*log(2*M_PI); //1.0/sqrt(2*M_PI);
	// for (int i = 0; i < x.size(); ++i)
	// {
	// 	double d = (x[i] - means[i])/errors[i];
	// 	v+= log_invSqrt - log(errors[i]) -0.5*d*d;
	// }
	// return log(peak) + v;

	double v = 0;
	// double logInv = -log(M_PI);
	for (int i =0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];
		v -= (log(1 + d*d) + log(errors[i]*M_PI)) ;
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
		// return -1.0/(errors[i] * errors[i]);
		double e = errors[i];
		double d = x[i] - means[i];
		double num = -2*(e - d) * (e+d);
		double denom = e*e + d*d;

		return num/(denom*denom);
	}
}

void dlogfunc(std::vector<double> & grad, const std::vector<double> & x)
{
	for (int i = 0; i < x.size(); ++i)
	{
		double d = (x[i] - means[i])/errors[i];
		// grad[i] = -1.0/(errors[i] * errors[i]);
		grad[i] = -2*d/(errors[i]* (1 + d*d));
	}
}
