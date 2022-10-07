#pragma once
#include <vector>
#include "targetFunction.h"
#include "funcs.h"


double test_MCI(int res, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<double> x(dim);
	double running = 0;
	// int useful =  0;
	for (int i = 0; i < res; ++i)
	{
		randomFill(x,lowerBound,upperBound);
		running += func(x);
	}

	double Volume = pow(upperBound - lowerBound,dim);
	return log(running * Volume / res);
}

double vegas_MCI(int res, int depth, int bins, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<std::vector<double>> hists(dim,std::vector<double>(bins,1.0/bins)); 
	std::vector<std::vector<double>> edges(dim,std::vector<double>(bins+1,1.0/bins)); 

	generateAccumulator(edges,hists,lowerBound,upperBound);
	double delta = (upperBound - lowerBound)/bins;
	double logDelta = log(delta);
	int reservedRes = res/2;
	int resPerDepth = std::max((res - reservedRes)/depth,5*bins);
	// std::cout << "Res per depth = " << resPerDepth << std::endl;
	std::vector<double> x(dim);
	std::vector<int> x_id(dim);
	double log_s = 0;

	double mem = 0.3;

	
	// std::vector<double> baseline = JSL::Vector::linspace(lowerBound,upperBound,bins);

	
	for (int i = 0; i < depth; ++i)
	{
		// std::fill(hists.begin(),hists.end(),std::vector<double>(bins,0.0));
		std::vector<std::vector<double>> new_hists(dim,std::vector<double>(bins,0)); 
		std::vector<std::vector<int>> hit_hist(dim,std::vector<int>(bins,0)); 
		for (int j = 0; j < resPerDepth; ++j)
		{
			// std::cout << "Generating a new position at " << j << ", bounded by " << lowerBound << "  " << upperBound << std::endl;
			// double w = fillFromDists(x,x_id,hists,lowerBound,upperBound);
			double w = fillFromDeltas(x,x_id,edges,lowerBound,upperBound);
			// std::cout << "\tres = " << j << "  " << JSL::Vector(x_id) << "  " << "| " << JSL::Vector(means) << "  weight = " << w << std::endl;
			double f = func(x);
			for (int k = 0; k < dim; ++k)
			{
				// if (hit_hist[k][x_id[k]])
				// {
				// 	new_hists[k][x_id[k]]=ale(new_hists[k][x_id[k]],f - w);
				// }
				// else
				// {
				// 	new_hists[k][x_id[k]]=f - w;
					
				// }
				new_hists[k][x_id[k]] += f/(exp(w)+1e-300);
				++hit_hist[k][x_id[k]];
			}
		}

		for (int j = 0; j < dim; ++j)
		{
			double run = 0;
			// std::cout << j << "\t" << JSL::Vector(hit_hist[j])<< std::endl;
			
			for (int b = 0; b < bins; ++b)
			{
				
				if (!hit_hist[j][b])
				{
					new_hists[j][b] = mem*hists[j][b] + 1e-200 ;
				}
				else
				{
					new_hists[j][b] /= hit_hist[j][b];
				}
			}
			if (dim > 30)
			{
				smooth(new_hists[j]);
			}
			for (int b = 0; b < bins; ++b)
			{
				// if (b == 0)
				// {
				// 	run = new_hists[j][b];
				// }
				// else
				// {
				// 	run = ale(run,new_hists[j][b]);
				// }
				run += new_hists[j][b];
			}
			for (int b = 0; b < bins; ++b)
			{
				hists[j][b] = mem * hists[j][b] + (1.0 - mem) * new_hists[j][b] /run;
				
			}
			// std::cout << "\t" << JSL::Vector(hists[j]) << std::endl;
		}
		// std::cout << "Depth = " <<i << " hist = " << JSL::Vector(hists[0]) << std::endl;
		generateAccumulator(edges,hists,lowerBound,upperBound);
		// if ((i+1)%1000 == 0 || i == depth-1)
		// {
		// 	JSL::gnuplot gp;
		// 	namespace lp = JSL::LineProperties;
		// 	for (int z = 0; z < dim; ++z)
		// 	{
		// 		// std::cout << " z = " << JSL::Vector(hists[z]) << std::endl;
		// 		gp.Plot(baseline,hists[z],lp::Legend("Dim" + std::to_string(z) + " at depth " + std::to_string(i)));
		// 	}
	
		// 	gp.SetLegend(true);
		// 	// gp.SetXRange(-3,3);
		// 	gp.SetYLog(true);
		// 	gp.Show();
		// }
	}
	
	
	

	double run = 0;
	for (int i = 0; i < reservedRes; ++i)
	{
		// double w = fillFromDists(x,x_id,hists,lowerBound,upperBound);
		double w = fillFromDeltas(x,x_id,edges,lowerBound,upperBound);
		{
			double f = func(x) / (exp(w)+1e-200);
			// if (i == 0)
			// {
			// 	run = f;
			// }
			// else
			// {
			// 	run = ale(run,f);
			// }
			run += f;
		}
	}
	return log(run /reservedRes);
}
