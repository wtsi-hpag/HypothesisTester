#pragma once
#include <random>

//very basic object. I could use a vector or something, but this keeps it obvious what is what
struct coordinate
{
	double x;
	double y;

};


double randomval(double lower, double upper)
{
	return lower + (double)rand()/RAND_MAX * (upper - lower);
}

std::vector<coordinate> generateData(int order, int nData)
{
	std::vector<coordinate> output(nData);

	//generate the coefficients of the polynomial
	std::vector<double> coefficient(order+1,0);
	for (int i = 0; i <= order; ++i)
	{
		coefficient[i] = randomval(-1,1);
	}

	double bounds = 3;
	double sigma = 1;
	std::random_device rd{};
    std::mt19937 gen{rd()};
	std::normal_distribution<> noise(0,sigma);
	
	//generate some random data points
	for (int i = 0; i < nData; ++i)
	{
		double x = randomval(-bounds,bounds);
		double y = 0;
		for (int p = 0; p <=order; ++p)
		{
			y += coefficient[p] * pow(x,p);
		}
		output[i].x = x;
		output[i].y = y + noise(gen); //add some noise in else the task is either trivial or impossble (dpeending on floating point accuracy)
	}
	return output;
} 