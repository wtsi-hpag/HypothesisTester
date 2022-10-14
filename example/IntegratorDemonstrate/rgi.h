#pragma once
#include <vector>
#include "targetFunction.h"
#include "funcs.h"

double test_RGI(int res, double lowerBound, double upperBound)
{
	int dim = means.size();

	int resPerDim = pow((double)res,1.0/dim);
	if (resPerDim < 2)
	{
		return log(0);
	}
	double delta = 1.0/(resPerDim-1);
	std::vector<double> x(dim);
	std::vector<int> dimProgress(dim,0);
	// std::vector<int> deltas(dim,delta);
	std::vector<int> progMax(dim,resPerDim);

	int total = pow(resPerDim,dim);
	
	// std::cout <<"Generating " << JSL::Vector(progMax) << std::endl;
	int bin;
	while (total < res)
	{
		bin = rand() % dim;
		progMax[bin] += 1;
		
		total = 1;
		for (int i = 0; i < dim; ++i)
		{
			total *= progMax[i];
		}
	}
	if (total > res)
	{
		progMax[bin]--;
		total *= (double)progMax[bin]/(progMax[bin] + 1);
	}

	// std::cout << "Without mods, would generate " <<pow(resPerDim,dim) << " now generating " << total << std::endl;

	// std::cout << JSL::Vector(progMax) << std::endl;
	double run = 0;
	int it = 0;
	bool integrationContinues = true;
	while (integrationContinues)
	{
		for (int i = 0; i < dim; ++i)
		{
			double xx = lowerBound + (double)dimProgress[i]/(progMax[i]-1) * (upperBound - lowerBound);
			x[i] = xx;
		}

		run += func(x);
		++it;

		//increment the vector
		int q = 0;
		++dimProgress[q];
		while (dimProgress[q] >= progMax[q])
		{
			if (q < dim - 1)
			{
				dimProgress[q] = 0;
				++dimProgress[q+1];
				++q;
			}
			else
			{
				integrationContinues = false;
				break;
			}
			// std::cout << "\t" << q << "  " << JSL::Vector(dimProgress) << std::endl;
		}
		// std::cout << JSL::Vector(dimProgress) << "   " << JSL::Vector(x) << std::endl;
	}

	return log(run/it)+ dim * log(upperBound - lowerBound);
}



double trapz_RGI(int res, double lowerBound, double upperBound)
{
	Status("tRGI",res);
	int dim = means.size();

	int resPerDim = pow((double)res,1.0/dim);
	if (resPerDim < 2)
	{
		return log(0);
	}
	double delta = 1.0/(resPerDim-1);
	std::vector<double> x(dim);
	std::vector<int> dimProgress(dim,0);
	// std::vector<int> deltas(dim,delta);
	std::vector<int> progMax(dim,resPerDim);

	int total = pow(resPerDim,dim);
	
	// std::cout <<"Generating " << JSL::Vector(progMax) << std::endl;
	int bin;
	while (total < res)
	{
		bin = rand() % dim;
		progMax[bin] += 1;
		
		total = 1;
		for (int i = 0; i < dim; ++i)
		{
			total *= progMax[i];
		}
	}
	if (total > res)
	{
		progMax[bin]--;
		total *= (double)progMax[bin]/(progMax[bin] + 1);
	}

	// std::cout << "Without mods, would generate " <<pow(resPerDim,dim) << " now generating " << total << std::endl;

	// std::cout << JSL::Vector(progMax) << std::endl;
	double run = 0;
	int it = 0;
	bool integrationContinues = true;
	while (integrationContinues)
	{
		double w = 1;
		for (int i = 0; i < dim; ++i)
		{
			double xx = lowerBound + (double)dimProgress[i]/(progMax[i]-1) * (upperBound - lowerBound);
			if (dimProgress[i] == 0 || dimProgress[i] == (progMax[i] -1))
			{
				w/=2;
			}
			x[i] = xx;
		}
		// std::cout << w << "   " << JSL::Vector(dimProgress) << std::endl;
		run += w * func(x);
		++it;

		//increment the vector
		int q = 0;
		++dimProgress[q];
		while (dimProgress[q] >= progMax[q])
		{
			if (q < dim - 1)
			{
				dimProgress[q] = 0;
				++dimProgress[q+1];
				++q;
			}
			else
			{
				integrationContinues = false;
				break;
			}
			// std::cout << "\t" << q << "  " << JSL::Vector(dimProgress) << std::endl;
		}
		// std::cout << JSL::Vector(dimProgress) << "   " << JSL::Vector(x) << std::endl;
	}
	Complete();
	return log(run/it)+ dim * log(upperBound - lowerBound);
}



double test_LRGI(int res, double lowerBound, double upperBound)
{
	int dim = means.size();

	int resPerDim = pow((double)res,1.0/dim);
	if (resPerDim < 2)
	{
		return log(0);
	}
	double delta = 1.0/(resPerDim-1);
	std::vector<double> x(dim);
	std::vector<int> dimProgress(dim,0);
	// std::vector<int> deltas(dim,delta);
	std::vector<int> progMax(dim,resPerDim);

	int total = pow(resPerDim,dim);
	
	// std::cout <<"Generating " << JSL::Vector(progMax) << std::endl;
	int bin;
	while (total < res)
	{
		bin = rand() % dim;
		progMax[bin] += 1;
		
		total = 1;
		for (int i = 0; i < dim; ++i)
		{
			total *= progMax[i];
		}
	}
	if (total > res)
	{
		progMax[bin]--;
		total *= (double)progMax[bin]/(progMax[bin] + 1);
	}

	// std::cout << "Without mods, would generate " <<pow(resPerDim,dim) << " now generating " << total << std::endl;

	// std::cout << JSL::Vector(progMax) << std::endl;
	double run = -1e300;
	int it = 0;
	bool integrationContinues = true;
	while (integrationContinues)
	{
		for (int i = 0; i < dim; ++i)
		{
			double xx = lowerBound + (double)dimProgress[i]/(progMax[i]-1) * (upperBound - lowerBound);
			x[i] = xx;
		}

		run = ale(logfunc(x),run);
		++it;

		//increment the vector
		int q = 0;
		++dimProgress[q];
		while (dimProgress[q] >= progMax[q])
		{
			if (q < dim - 1)
			{
				dimProgress[q] = 0;
				++dimProgress[q+1];
				++q;
			}
			else
			{
				integrationContinues = false;
				break;
			}
			// std::cout << "\t" << q << "  " << JSL::Vector(dimProgress) << std::endl;
		}
		// std::cout << JSL::Vector(dimProgress) << "   " << JSL::Vector(x) << std::endl;
	}
	// double Volume = pow(upperBound - lowerBound,dim);
	// std::cout << ""
	// std::cout << "Wanted " << res << " got " << it << std::endl;
	return run - log(it) + dim * log(upperBound - lowerBound);
}