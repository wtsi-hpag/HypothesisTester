#pragma once
#include <vector>
#include <iostream>
#include <math.h>
#include <string>
#include "matrix.h"
double explogadd(double loga, double logb)
{
	//An efficient way to compute log(a+b) given log(a) and log(b), whilst ensuring as few exp-overflows as possible
	if (loga > logb)
	{
		return loga + log(1 + exp(logb - loga));
	}
	else
	{
		return logb + log(1 + exp(loga - logb));
	}
};


enum integrator {MCI,LMCI,GAI};
template<class DataClass>
class Hypothesis
{
	public:
		Hypothesis(int dimension) : Dimension(dimension)
		{
			paramLowerBound.resize(Dimension,-1000);
			paramUpperBound.resize(Dimension,1000);
			uniformPriorComputed= false;
			Mode = LMCI;
			Identifier = "Undifferentiated Hypothesis of dimension " + std::to_string(dimension);
		}
		Hypothesis(int dimension, integrator mode) : Dimension(dimension)
		{
			paramLowerBound.resize(Dimension,-1000);
			paramUpperBound.resize(Dimension,1000);
			uniformPriorComputed= false;
			Mode = mode;
			Identifier = "Undifferentiated Hypothesis of dimension " + std::to_string(dimension);
		}
		void SetLowerBound(const std::vector<double> lowerBound)
		{
			assert(lowerBound.size() == Dimension);
			paramLowerBound = lowerBound;
		}
		void SetUpperBound(const std::vector<double> bound)
		{
			assert(bound.size() == Dimension);
			paramUpperBound = bound;
		}
		
		virtual std::vector<double> FindMaximum(const std::vector<DataClass> & Data, std::vector<double> & initialGuess, int N)
		{
			double alpha = 1e-3;
			double alphaMax = 1e2;
			int accelSteps = N/20;
			double scaling = pow(alphaMax/alpha,1.0/accelSteps);
			// int N = 10000;
			std::vector<double> grad(Dimension);
			std::vector<double> pos = initialGuess;
			std::vector<double> m(Dimension,0);
			std::vector<double> v(Dimension,0);
			double b1_init = 0.3;
			double b2_init = 0.6;
			double b1_true = 0.7;
			double b2_true = 0.999;
			double eps = 1e-8;
			double prevValue = 0;
			double nablaThreshold = sqrt(Dimension) * 1e-10;
			double absGrad;
			double minVal = 0;
			std::vector<double> is;
			std::vector<double> absGrads;
			std::vector<double> a0s;
			std::vector<double> vals;
			std::vector<double> sigs;
			std::vector<double> previousPosition;
			bool haltRecompute = false;
			double functionVal = LogProbability(Data,pos) + LogPrior(pos);
			int burnInSteps = 50;
			int t = 0;
			double b1 = b1_init;
			double b2 = b2_init;
			for (int i = 0; i < N+burnInSteps; ++i)
			{
				// std::cout <<"optimisation step " << i << "/" << N + burnInSteps<< std::endl;
				if (i == burnInSteps)
				{
					std::fill(m.begin(),m.end(),0.0);
					std::fill(v.begin(),v.end(),0.0);
					t = 0;
					b1 = b1_true;
					b2 = b2_true;
					alpha = alpha/10;
				}
				++t;

				if (haltRecompute == false)
				{
					LogProbabilityGradient(grad,Data,pos);
				
					absGrad = 0;
					for (int i = 0; i < grad.size(); ++i)
					{
						absGrad += grad[i]*grad[i];
					}
					absGrad = sqrt(absGrad);
				}
				// std::cout << i << "/" << N + burnInSteps << "  " << JSL::Vector(pos) << "   " << JSL::Vector(grad) << "  " << absGrad << "  " << alpha << "  " << functionVal << std::endl;
				if (i >= 0)
				{
					is.push_back(i);
					absGrads.push_back(absGrad);
					a0s.push_back(pos[Dimension-2]);
					vals.push_back(functionVal);
					sigs.push_back(pos[Dimension-1]);
				}
				if (i == 0 || absGrad < minVal)
				{
					minVal = absGrad;
				}
				if (absGrad < nablaThreshold && i >= burnInSteps)
				{
					break;
				}
				else
				{
				// grad[Dimension-1] = 0;
					for (int j = 0; j < Dimension; ++j)
					{
						m[j] = b1 * m[j] + (1.0 - b1)*grad[j];
						v[j] = b2* v[j] + (1.0 - b2) * grad[j]*grad[j];

						double correctedM = m[j]/(1.0 - pow(b1,t));
						double correctedV = v[j]/(1.0 - pow(b2,t));


						pos[j] += alpha * correctedM/(sqrt(correctedV) + eps);
						pos[j] = std::min(paramUpperBound[j],std::max(paramLowerBound[j],pos[j]));
					}

					double pp = LogProbability(Data,pos) + LogPrior(pos);

					
					if (i == 0 || pp >= prevValue)
					{
						alpha = alpha * scaling;
						haltRecompute = false;
						// previousPosition = pos;
					}
					else
					{
						// haltRecompute = true;
						// pos = previousPosition;
						alpha/=pow(scaling,5);
						// pp = prevValue;
						
					}
					prevValue = pp;
					functionVal = pp;
				}
			}
			return pos;
		}
	
		std::vector<double> paramLowerBound;
		std::vector<double> paramUpperBound;
		std::string Identifier;

		std::vector<double> paramMidPoint()
		{
			std::vector<double> mid(Dimension);
			for (int i = 0; i < Dimension; ++i)
			{
				mid[i] = 0.5 * (paramLowerBound[i] + paramUpperBound[i]);
			}
			return mid;
		}

		virtual double Score(const std::vector<DataClass> Data,int resolution)
		{
			PrecomputeData(Data);
			if (Dimension == 0)
			{
				return LogProbability(Data,{});
			}
			switch (Mode)
			{
				case LMCI:
				{
					double s = LMCI_Implementation(Data,resolution);
					return s;
				}
				case MCI:
				{
					return MCI_Implementation(Data,resolution);
				}
				case GAI:
				{
					return GAI_Implementation(Data,resolution);
				}
				default:
				{
					std::cout << "ERROR: No support for chosen integrator, using MonteCarlo instead" << std::endl;
					return MCI_Implementation(Data,resolution);
				}
			}
		}

		int Verbosity = 0;
		double Tolerance = 0.00000;
		int Dimension;
	private:
			
		bool uniformPriorComputed;
		double UniformPriorValue;
		int iterations = 0;
		integrator Mode;
		double IntegrationVolume()
		{
			double V = 1;
			double logV = 0;
			for (int i = 0;i < Dimension; ++i)
			{
				
				V *= (paramUpperBound[i] - paramLowerBound[i]);
				logV += log(paramUpperBound[i] - paramLowerBound[i]);
			}
			UniformPriorValue = -logV;
			uniformPriorComputed = true;
			return V;
		}

		virtual void PrecomputeData(const std::vector<DataClass> & Data)
		{
		
		}



		virtual double LogProbability(const std::vector<DataClass> & Data, const std::vector<double> & params)
		{	
			std::cout <<"Model Probability Called from unspecialised superclass - you must provide an override for the definition of LogProbability" << std::endl;
			exit(1);
			return 0;
		};
		
		
		
		virtual double LogPrior(const std::vector<double> & params)
		{
			if (uniformPriorComputed)
			{
				return UniformPriorValue;
			}
			else
			{
				uniformPriorComputed = true;
				IntegrationVolume();
				// UniformPriorValue = -1 * log(IntegrationVolume());
				return UniformPriorValue;
			}
		};
	
		virtual void LogProbabilityGradient(std::vector<double> & g, const std::vector<DataClass> & Data, std::vector<double> & params)
		{
			double v = LogProbability(Data,params);
			for (int dim = 0; dim < Dimension; ++dim)
			{
				double orig = params[dim];
				double delta = std::max(1e-8,abs(params[dim]*1e-2));

				
				params[dim] += delta;
				double vPlus = LogProbability(Data,params);
				params[dim] = orig;

				g[dim] = (vPlus - v)/delta;
			}
		}


		virtual squareMatrix ProbabilityHessian(const std::vector<DataClass> & Data, std::vector<double> & params)
		{
			std::cout << "No default Hessian implemented yet" << std::endl;
			squareMatrix m(Dimension);
			return m;
		}
			double LMCI_Implementation(const std::vector<DataClass> & Data, int resolution)
		{
			double s = 0;
			std::vector<double> params(Dimension);
			double V = IntegrationVolume();
			for (int i = 0; i < resolution; ++i)
			{
				double rSq = 0;
				for (int j = 0; j < Dimension; ++j)
				{	
					double r = (double)rand()/RAND_MAX;
					params[j] = paramLowerBound[j] + r * (paramUpperBound[j] - paramLowerBound[j]);
				}
				double pp = LogProbability(Data,params);
				double logp = LogPrior(params) + pp;
				
				if (i == 0)
				{
					s = logp;
				}
				else
				{
					s = explogadd(s,logp);
				}
				// std::cout << "Integral step " << i << ":  " << params[0] << "  " << logp << "  " << s << std::endl;
			}
			
			double score = -1*UniformPriorValue - log(resolution) + s;
			return score;
		}
		double MCI_Implementation(const std::vector<DataClass> & Data, int resolution)
		{
			double s = 0;
			std::vector<double> params(Dimension);
			double V = IntegrationVolume();
			if (Verbosity >=1)
			{
				std::cout << "\t\tNow beginning " << resolution << " samples of the function" << std::endl;
			}
			for (int i = 0; i < resolution; ++i)
			{
				double rSq = 0;
				for (int j = 0; j < Dimension; ++j)
				{	
					double r = (double)rand()/RAND_MAX;
					params[j] = paramLowerBound[j] + r * (paramUpperBound[j] - paramLowerBound[j]);
				}
				double pp = exp(LogProbability(Data,params));
				double p = exp(LogPrior(params)) * pp;

				s += p;
			}
		
			double score = log(V/resolution) + log(s);
		
			return score;
		}

		double GAI_Implementation(const std::vector<DataClass> & Data, int resolution)
		{
			auto start = paramMidPoint();
			auto optimum = FindMaximum(Data,start,resolution);
			double optimalValue = LogProbability(Data,optimum);
			if (Verbosity >=1)
			{
				std::cout << "\t\tOptimum value = " << optimalValue << "\n\t\tNow computing det(H)" << std::endl;
			}
			squareMatrix s = ProbabilityHessian(Data,optimum);
			double d= s.log_LU_Determinant();
			int inferredDimension = s.Dimension; //allows individual hypotheses to lower the dimension of the Hessian, i.e. if a bunch of variables are superfluous and don't enter into the Gaussian form
			if (Verbosity >=1)
			{
				std::cout << "\t\t" << inferredDimension << "-dimension Hessian computed with determinant " << d << std::endl;
			}
			double v = optimalValue +(double)inferredDimension/2 * log(2*M_PI) -0.5 * d;
			
			return v;
			
		}
};
