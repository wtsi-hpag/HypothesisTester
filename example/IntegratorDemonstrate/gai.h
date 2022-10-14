#pragma once
#include <vector>
#include "targetFunction.h"
#include "../../HypothesisTester.h"
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

void manualGrad(std::vector<double> & grad,std::vector<double> x)
{
	double delta_base = 1e-3;
	for (int i = 0; i < x.size(); ++i)
	{
		double orig = x[i];
		double delta = std::max(1e-10,delta_base * orig);

		x[i] += delta;
		double up = logfunc(x);
		x[i] = orig - delta;
		double down = logfunc(x);
		x[i] = orig;

		grad[i] = (up - down)/(2*delta);
	}
}
squareMatrix manualHessian(std::vector<double> x)
{
	int d = x.size();
	squareMatrix out(d);

	double base = logfunc(x);
	double deltaBase = 1e-2;
	for (int i = 0; i < d; ++i)
	{
		double orig_x = x[i];
		double deltaX = std::max(1e-5,deltaBase * orig_x);

		x[i] = orig_x + deltaX;
		double upX = logfunc(x);
		x[i] = orig_x - deltaX;
		double downX = logfunc(x);
		x[i] = orig_x;

		double fxx = (upX + downX - 2*base)/(deltaX * deltaX);
		out(i,i) = fxx;

		for (int j = i+1; j < d; ++j)
		{
			double orig_y = x[j];
			double deltaY = std::max(1e-5,deltaBase * orig_y);

			x[i] = orig_x + deltaX;
			x[j] = orig_y + deltaY;
			double upXupY = logfunc(x);
			x[i] = orig_x;
			double upY = logfunc(x);
			x[j] = orig_y - deltaY;
			double downY = logfunc(x);
			x[i] = orig_x - deltaX;
			double downXdownY = logfunc(x);

			double fxy = (upXupY - upX - upY + 2*base - downX - downY + downXdownY)/(deltaX * deltaY);
			out(i,j) = fxy;
			out(j,i) = fxy;
			x[i] = orig_x;
			x[j] = orig_y;
		}

	}
	return out;
}

double gaiSolve(const std::vector<double> pos)
{
	
	double f = logfunc(pos);

	squareMatrix H = Hessian(pos);
	
	// std::cout << "f = " << f << std::endl;
	// std::cout << "H = " << H.log_LU_Determinant() << "" << std::endl;
	// std::cout << "\nVal = " << f + exact.size()* log(2*M_PI)/2 - 0.5 * H.log_LU_Determinant() << "\n\n";
	return f + pos.size()* log(2*M_PI)/2 - 0.5 * H.log_LU_Determinant();
}

double exact_GAI(int res, double lowerBound, double upperBound)
{
	Status("GAI-EXACT",0);
	double v = gaiSolve(means);
	Complete();
	return v;
}


double test_GAI(int res, double lowerBound, double upperBound)
{
	Status("GAI-FUNCTION",res);
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
		if (i % 5 == 0)
		{
			for (int j = 0; j < d; ++j)
			{
				guess[j] = std::min(upperBound,std::max(lowerBound,guess[j]));
			}
		}
		// std::cout << "Corrected to: = " << guess << std::endl;
	}

	// std::cout << "Guess: " << guess << "\nActual: " << JSL::Vector(means) << std::endl;
	// std::cout << "Momentum: " << momentum << "\nGrad: " << JSL::Vector(grad) << std::endl;
	// std::cout << "Final alpha = " << alpha << std::endl;

	double v =  gaiSolve(guess);
	Complete();
	return v;
}

double basic_GAI(int res, double lowerBound, double upperBound)
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
		// dlogfunc(grad,guess);
		manualGrad(grad,guess);
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
		if (i%5 == 0)
		{
			for (int j = 0; j < d; ++j)
			{
				guess[j] = std::min(upperBound,std::max(lowerBound,guess[j]));
			}
		}
	}

	return logfunc(guess) + guess.Size()* log(2*M_PI)/2 - 0.5 * manualHessian(guess).log_LU_Determinant();
}