#pragma once
#include <vector>
#include <string>
#include "Hypothesis.h"
#include <thread>
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
		ModelTester()
		{
			Nthreads = 1;
		}

		ModelTester(int n)
		{
			Nthreads = n;
		}

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


			Scores = std::vector(N,0.0);
			Threads.resize(Nthreads-1);
			Names.resize(N);

			std::vector<int> count(Nthreads,0);
			int allocated = 0;
			int idx = 0;
			while (allocated < N)
			{
				++count[idx % Nthreads];
				++allocated;
				++idx;
			}

			int t = 0;
			for (int i = 0; i < Nthreads -1; ++i)
			{
				Threads[i] = std::thread(&ModelTester::ChunkedScoreLauncher,this,Data,resolution,t,count[i]);
				t+=count[i];
			}
			// Threads[Nthreads-1]= std::thread(&ModelTester::ChunkedScoreLauncher,this,Data,resolution,t,count[Nthreads-1]);
			ChunkedScoreLauncher(Data,resolution,t,count[Nthreads-1]);

			int joined = 1;
			while (joined < Nthreads)
			{
				for (int n = 0; n < Nthreads-1; ++n)
				{
					if (Threads[n].joinable())
					{
						Threads[n].join();
						++joined;
					}
				}
			}

			double bestScore;
			int bestHyp;
			for (int i = 0; i <N; ++i)
			{
				if (i==0 || Scores[i] > bestScore)
				{
					bestScore = Scores[i];
					bestHyp = i;
				}
			}


			for (int i = 0; i < Suppositions.size(); ++i)
			{
				Scores[i] -= bestScore;
			}
			InferenceResults output;
			output.Models = Names;
			output.Scores = Scores;
			output.BestModel = Names[bestHyp];
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
		bool PrintMessages = false;
	private:
		std::vector<std::unique_ptr<Hypothesis<DataClass>>> Suppositions; 
		std::vector<double> Scores;
		std::vector<std::thread> Threads;
		std::vector<std::string> Names;
		int Nthreads;

		void ChunkedScoreLauncher(const std::vector<DataClass> & Data, int resolution,int start, int size)
		{
			int end = std::min(start+size,(int)Suppositions.size());
			for (int i = start; i < end; ++i)
			{
				if (PrintMessages)
				{
					std::cout << "Beginning Test on " << Suppositions[i]->Identifier << std::endl;
				}
				double S = Suppositions[i]->Score(Data,resolution);
				Scores[i] = S;
				Names[i] =Suppositions[i]->Identifier;
			}
		}
};