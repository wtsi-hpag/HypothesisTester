#pragma once 
#include "../Hypothesis.h"
#include "coordinate.h"


/*
	This is an example of how to create suitable entities for the HypothesisTester engine to act upon.
	All classes should be child classes of Hypothesis<data> class, and -- at the absolute bare minimum, should declare a LogProbability() function as below.
	If a suitabe LogProbability is not defined, the code will forcibly exit, as the hypothesis is therefore meaningless. 
*/
class PolynomialHypothesis : public Hypothesis<coordinate>
{
	public:
		PolynomialHypothesis(int order) : Hypothesis<coordinate>(order+2)
		{
			Order = order;
			//2 additional parameters come from the zeroth order parameter, and from the error which we are also inferring
			std::vector<double> lower(order+2,-2);
			std::vector<double> upper(order+2, 2);
			lower[order+1] = -5; //error has different bounds than the coefficients!
			upper[order+1]= 1; 

			SetLowerBound(lower);
			SetUpperBound(upper);
			Identifier = "Polynomial of Order " + std::to_string(order); //gives it a nice name to be printed out at the end
		}


		double LogProbability(const std::vector<coordinate> & Data, const std::vector<double> & params)
		{	
			double logSigma = params[Order+1];
			double sigma = exp(logSigma);

			double L = 0;

			for (int i = 0; i < Data.size(); ++i)
			{
				double predictedY = 0;
				for (int p = 0; p <= Order; ++p)
				{
					predictedY += params[p] * pow(Data[i].x,p);
				}
				double d = (predictedY - Data[i].y)/sigma;
				L -= (logSigma + 0.5 * d*d);
			}

			return L;
		};

	private:
		int Order;




};