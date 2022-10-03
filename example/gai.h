#pragma once
#include <vector>
#include "targetFunction.h"
#include "../HypothesisTester.h"
#include "JSL.h"

squareMatrix Hessian(const std::vector<double> & x)
{
	int d = x.size();
	squareMatrix out(d);
	for (int i = 0; i < d; ++i)
	{
		for (int j = i; j < d; ++j)
		{
			out(i,j) = ddlogfunc(x,i,j);
			out(j,i) = out(i,j);
		}

	}
	return out;
}




double cheat_GAI(const std::vector<double> exact)
{
	double f = logfunc(exact);

	squareMatrix H = Hessian(exact);
	
	// std::cout << "f = " << f << std::endl;
	// std::cout << "H = " << H.log_LU_Determinant() << "" << std::endl;
	// std::cout << "\nVal = " << f + exact.size()* log(2*M_PI)/2 - 0.5 * H.log_LU_Determinant() << "\n\n";
	return f + exact.size()* log(2*M_PI)/2 - 0.5 * H.log_LU_Determinant();
}

double test_GAI(int res, double lowerBound, double upperBound)
{
	int d = means.size();
	JSL::Vector guess(d);
	guess += (upperBound + lowerBound)/2;
	JSL::Vector momentum(d);
	double mem = 0.6;
	std::vector<double> grad(d);
	double alpha = 1;
	double prevF;
	for (int i = 0; i < res; ++i)
	{
		dlogfunc(grad,guess);

		momentum = mem*momentum + (1.0 - mem)* JSL::Vector(grad);
		
		guess += alpha * momentum;

		double newF = logfunc(guess);
		if (i > 0)
		{
			if (newF > prevF)
			{
				alpha *= 1.05;
			}
			else
			{
				alpha /= 10;
			}
		}
		prevF = newF;
		// std::cout << "Target = " << JSL::Vector(means) << std::endl;
		// std::cout << "Guess = " << guess << std::endl;
		for (int i = 0; i < d; ++i)
		{
			guess[i] = std::min(upperBound,std::max(lowerBound,guess[i]));
		}
		// std::cout << "Corrected to: = " << guess << std::endl;
	}

	// std::cout << "Guess: " << guess << "\nActual: " << JSL::Vector(means) << std::endl;
	// std::cout << "Momentum: " << momentum << "\nGrad: " << JSL::Vector(grad) << std::endl;
	// std::cout << "Final alpha = " << alpha << std::endl;
	return cheat_GAI(guess);
}