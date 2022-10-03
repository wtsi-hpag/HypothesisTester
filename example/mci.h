#pragma once
#include <vector>
#include "targetFunction.h"
#include "funcs.h"


double test_MCI(int res, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<double> x(dim);
	double running = 0;
	// int useful =  0;
	for (int i = 0; i < res; ++i)
	{
		randomFill(x,lowerBound,upperBound);
		running += func(x);
	}

	double Volume = pow(upperBound - lowerBound,dim);
	return running * Volume / res;
}
