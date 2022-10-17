#pragma once
#include "funcs.h"
#include "targetFunction.h"
#include "JSL.h"

double test_LMCI(int res, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<double> x(dim);
	double running = 0;
	for (int i = 0; i < res; ++i)
	{
		randomFill(x,lowerBound,upperBound);
		double t = logfunc(x);
		if (i == 0)
		{
			running = t;
		}
		else
		{
			running = ale(running,t);
		}
	}

	double logVolume = dim *log(upperBound - lowerBound);
	return running + logVolume - log(res);
}



double vegas_LMCI(int res, int depth, int bins, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<std::vector<double>> hists(dim,std::vector<double>(bins,1.0/bins)); 
	std::vector<std::vector<double>> edges(dim,std::vector<double>(bins+1,1.0/bins)); 

	generateAccumulator(edges,hists,lowerBound,upperBound);

	int reservedRes = res/2;
	int resPerDepth = std::max((res - reservedRes)/depth,10);
	// std::cout << "Res per depth = " << resPerDepth << std::endl;
	std::vector<double> x(dim);
	std::vector<int> x_id(dim);

	double mem = 0.3;

	
	std::vector<std::vector<double>> new_hists(dim,std::vector<double>(bins,-999999999)); 
	std::vector<std::vector<int>> hit_hist(dim,std::vector<int>(bins,0)); 
	for (int i = 0; i < depth; ++i)
	{
		for (int q = 0; q < dim; ++q)
		{
			std::fill(new_hists[q].begin(),new_hists[q].end(),-1e200);
			std::fill(hit_hist[q].begin(),hit_hist[q].end(),0);
		}

		for (int j = 0; j < resPerDepth; ++j)
		{
			
			double w = fillFromDeltas(x,x_id,edges,lowerBound,upperBound);
			double f = logfunc(x);
			for (int k = 0; k < dim; ++k)
			{
				new_hists[k][x_id[k]]=ale(new_hists[k][x_id[k]],f - w);
				++hit_hist[k][x_id[k]];
			}
		}

		for (int j = 0; j < dim; ++j)
		{
			double run = 0;
			
			for (int b = 0; b < bins; ++b)
			{
				
				if (!hit_hist[j][b])
				{
					new_hists[j][b] = log(mem*hists[j][b] + 1e-200) ;
				}
				else
				{
					new_hists[j][b] -= log(hit_hist[j][b]);
				}
			}
			if (dim > 30)
			{
				smooth(new_hists[j]);
			}
			for (int b = 0; b < bins; ++b)
			{
				if (b == 0)
				{
					run = new_hists[j][b];
				}
				else
				{
					run = ale(run,new_hists[j][b]);
				}
			}
			for (int b = 0; b < bins; ++b)
			{
				hists[j][b] = mem * hists[j][b] + (1.0 - mem) * exp(new_hists[j][b] - run);
				
			}
		}
		generateAccumulator(edges,hists,lowerBound,upperBound);
	}

	double run = 0;
	for (int i = 0; i < reservedRes; ++i)
	{
		// double w = fillFromDists(x,x_id,hists,lowerBound,upperBound);
		double w = fillFromDeltas(x,x_id,edges,lowerBound,upperBound);
		{
			double f = logfunc(x) - w;
			if (i == 0)
			{
				run = f;
			}
			else
			{
				run = ale(run,f);
			}
		}
	}
	Complete();
	return run - log(reservedRes);
}
