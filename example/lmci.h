#pragma once
#include "funcs.h"
#include "targetFunction.h"
#include "JSL.h"

double test_LMCI(int res, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<double> x(dim);
	double running = 0;
	// int useful =  0;
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
		// running += func(x);
	}

	double logVolume = dim *log(upperBound - lowerBound);
	// std::cout << logVolume << std::endl;
	return running + logVolume - log(res);
}



double vegas_MCI(int res, int depth, int bins, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<std::vector<double>> hists(dim,std::vector<double>(bins,1.0/bins)); 
	double delta = (upperBound - lowerBound)/bins;
	double logDelta = log(delta);
	int reservedRes = res/2;
	int resPerDepth = std::max((res - reservedRes)/depth,5*bins);
	// std::cout << "Res per depth = " << resPerDepth << std::endl;
	std::vector<double> x(dim);
	std::vector<int> x_id(dim);
	double log_s = 0;

	double mem = 0.6;

	
	std::vector<double> baseline = JSL::Vector::linspace(lowerBound,upperBound,bins);

	
	for (int i = 0; i < depth; ++i)
	{
		// std::fill(hists.begin(),hists.end(),std::vector<double>(bins,0.0));
		std::vector<std::vector<double>> new_hists(dim,std::vector<double>(bins,-999999999)); 
		std::vector<std::vector<bool>> hit_hist(dim,std::vector<bool>(bins,false)); 
		for (int j = 0; j < resPerDepth; ++j)
		{
			// std::cout << "Generating a new position at " << j << ", bounded by " << lowerBound << "  " << upperBound << std::endl;
			double w = fillFromDists(x,x_id,hists,lowerBound,upperBound);

			// std::cout << "\tres = " << j << "  " << JSL::Vector(x_id) << "  " << "| " << JSL::Vector(means) << "  weight = " << w << std::endl;
			double f = logfunc(x);
			for (int k = 0; k < dim; ++k)
			{
				if (hit_hist[k][x_id[k]])
				{
					new_hists[k][x_id[k]]=ale(new_hists[k][x_id[k]],f - w);
				}
				else
				{
					new_hists[k][x_id[k]]=f - w;
					hit_hist[k][x_id[k]] = true;
				}
			}
		}

		for (int j = 0; j < dim; ++j)
		{
			double run = 0;
			for (int b = 0; b < bins; ++b)
			{

				if (!hit_hist[j][b])
				{
					// std::cout << "smoothing at " << j << "  " << b << std::endl;
					int low=-1;
					int high=-1;

					for (int l = b -1; l >= 0; --l)
					{
						if (hit_hist[j][l])
						{
							low = l;
							break;
						}
					}
					// std::cout << low << std::endl;
					for (int u = b + 1; u < hit_hist[j].size(); ++u)
					{
						if (hit_hist[j][u])
						{
							high = u;
							break;
						}
					}
					// std::cout <<high << std::endl;
					double bruch;
					double base;
					if (high >= 0 && low >= 0)
					{
						base = new_hists[j][low];
						bruch = (new_hists[j][high] - new_hists[j][low])/(high - low);
					}
					else
					{
						int active = std::max(high,low);
						base = new_hists[j][active];
						bruch = 0;
					}

					for (int bb = b; bb < std::min(high,(int)new_hists[j].size()); ++bb)
					{
						new_hists[j][bb] = base + bruch * (bb- low);
						// std::cout << "Interp at " << bb << std::endl;
						hit_hist[j][bb] = true;
					}
				}


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
		// std::cout << "Depth = " <<i << " hist = " << JSL::Vector(hists[0]) << std::endl;
		
		if (i%10 == 0 || i == depth-1)
		{
			JSL::gnuplot gp;
			namespace lp = JSL::LineProperties;
			for (int z = 0; z < dim; ++z)
			{
				// std::cout << " z = " << JSL::Vector(hists[z]) << std::endl;
				gp.Plot(baseline,hists[z],lp::Legend("Dim" + std::to_string(z) + " at depth " + std::to_string(i)));
			}
	
			gp.SetLegend(true);
			gp.SetXRange(-3,3);
			gp.Show();
		}
	}
	
	
	

	double run = 0;
	for (int i = 0; i < reservedRes; ++i)
	{
		double w = fillFromDists(x,x_id,hists,lowerBound,upperBound);
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
	return run - log(reservedRes);
}


double vegas2_MCI(int res, int depth, int bins, double lowerBound, double upperBound)
{
	int dim = means.size();

	std::vector<std::vector<double>> deltas(dim,std::vector<double>(bins,(upperBound - lowerBound)/bins)); 
	std::vector<std::vector<double>> accum(dim,std::vector<double>(bins+1,(upperBound - lowerBound)/bins)); 
	int reservedRes = res/2;
	int resPerDepth = std::max((res - reservedRes)/depth,5*bins);
	std::vector<double> x(dim);
	std::vector<int> x_id(bins);
	//burn in


	for (int d = 0; d < depth; ++d)
	{
		accumulate(deltas,accum,lowerBound,upperBound);
		std::vector<std::vector<double>> hist(dim,std::vector<double>(bins,-1e200)); 
		for (int t = 0; t < resPerDepth; ++t)
		{
			double w = fillFromDeltas(x,x_id,accum);
			double f = logfunc(x) - w;

			for (int i = 0; i < dim; ++i)
			{
				double & tgt = hist[i][x_id[i]];

				tgt = ale(tgt,f);
			}
		}
		for (int z = 0; z < dim; ++z)
			std::cout << "z = " << JSL::Vector(hist[z]) << std::endl;
		

		if (d%10 == 0 || d == depth-1)
		{
			JSL::gnuplot gp;
			namespace lp = JSL::LineProperties;
			std::vector<double> x(accum[0].size()*4);
			std::vector<double> y = x;
			for (int z = 0; z < dim; ++z)
			{
				// std::cout << " z = " << JSL::Vector(hists[z]) << std::endl;
				double acc = 0;
				for (int q = 0; q < accum[z].size()-1; ++q)
				{
					int v = 4*q;
					double w = 1.0/(bins * deltas[z][q]);
					x[v]=accum[z][q];
					x[v+1] = accum[z][q];
					x[v+2]=accum[z][q+1];
					x[v+3] = accum[z][q+1];
					y[v] = 0;
					y[v+1] = w;
					y[v+2] = w;
					y[v+3] = 0;
					acc += w * deltas[z][q];
				}
				std::cout << "Sum = " << acc << std::endl;
				// std::cout << "x = " << JSL::Vector(x) << " \ny= " << JSL::Vector(y) << std::endl;
				gp.Plot(x,y,lp::Legend("Dim" + std::to_string(z) + " at depth " + std::to_string(d)));
			}
	
			gp.SetLegend(true);
			// gp.SetXRange(-3,3);
			gp.Show();
		}
	}

	//actual integral

	accumulate(deltas,accum,lowerBound,upperBound);
	// std::cout << JSL::Vector(accum[0])<<std::endl;
	double run = 0;
	for (int i = 0; i < reservedRes; ++i)
	{
		double w = fillFromDeltas(x,x_id,accum);
		// std::cout << JSL::Vector(x) << " with weight " << w << std::endl;
		double f = logfunc(x) + w;
		if (i == 0)
		{
			run = f;
		}
		else
		{
			run = ale(run,f);
		}
		
	}
	return run - log(reservedRes);
}