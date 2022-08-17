#pragma once
#include <vector>
#include <string>
#include "Hypothesis.h"

struct InferenceResults
{
	std::vector<std::string> Models;
	std::vector<double> Scores;
	int BestModelIdx;
	std::string BestModel;
};


template<class DataClass>
class ModelTester
{
	public:
		ModelTester(){};

		template<class T>
		void AddHypothesis(T guess)
		{			
			std::unique_ptr<T> p = std::make_unique<T>(guess);
			Suppositions.push_back(std::move(p));
		}

		InferenceResults BeginTest(const std::vector<DataClass> & Data, int resolution)
		{
			int N = Suppositions.size();
			if (N < 2)
			{
				std::cout << "Cannot perform test with fewer than 2 hypotheses" << std::endl;
				exit(5);
			} 
			double bestScore = 0;
			int bestHyp = -1;
			std::vector<double> scores;
			std::vector<std::string> names;
			
			for (int i = 0; i < Suppositions.size(); ++i)
			{

				std::cout << "Testing " << Suppositions[i]->Identifier << std::endl;
				double S =Suppositions[i]->Score(Data,resolution);
				std::cout << "\tScored " << S << std::endl;
				if (i == 0 || (!std::isnan(S) && !std::isinf(S) && S > bestScore))
				{
					bestScore = S;
					bestHyp = i;
					std::cout << "\tAssigining best score" << std::endl;
				}
				std::cout << "\tBeginning model fit" << std::endl;
				scores.push_back(S);
				names.push_back(Suppositions[i]->Identifier);
				std::cout << "\tCompleted" << std::endl;
			}
			
			for (int i = 0; i < Suppositions.size(); ++i)
			{
				scores[i] -= bestScore;
			}
			InferenceResults output;
			output.Models = names;
			output.Scores = scores;
			output.BestModel = names[bestHyp];
			output.BestModelIdx = bestHyp;

			return output;
		}

		std::vector<double> FitModel(int modelIdx,const std::vector<DataClass> & Data)
		{
			
			auto start = Suppositions[modelIdx]->paramMidPoint();
			if (start.size() > 0)
			{
				return Suppositions[modelIdx]->FindMaximum(Data,start,100000);
			}
			else
			{
				return {};
			}
		}
	std::vector<std::unique_ptr<Hypothesis<DataClass>>> Suppositions; 
};