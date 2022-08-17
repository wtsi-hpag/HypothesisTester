#include "../HypothesisTester.h"
#include <iostream>
#include <vector>

#include "coordinate.h"
#include "PolynomialHypothesis.h"



/* 
	This is an example to demonstrate how the inference engine functions

	In coordinate.h we generate a nuber of (x,y) pairs according to a randomly generated polynomial of order N
	The goal is for the inference engine to infer the value of N from those (x,y) pairs

	This is achieved by PolynomialHypothesis objects -- defined in PolynomialHypothesis.h, each associated with a value of N.
*/



int main(int argc, char ** argv)
{

	unsigned int order = 4;
	unsigned int nData = 100;
	srand(time(NULL));
	
	//Allow basic command line interaction
	if (argc > 1)
	{
		order = std::stoi(argv[1]);
	}
	if (argc > 2)
	{
		nData = std::stoi(argv[2]);
	}
	
	auto data = generateData(order,nData);

	ModelTester<coordinate> T(true);

	//add one hypothesis at a time into the tester object
	for (int i = 0; i < order + 4; ++i)
	{
		T.AddHypothesis(PolynomialHypothesis(i));
	}

	//This is the only place where data can be introduced, to prevent statistical mishandling.
	auto results = T.BeginTest(data,100000);


	std::cout << "Best Fitting Model is: " << results.BestModel << std::endl;
	return 0;
}
