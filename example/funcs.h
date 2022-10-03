#pragma once
#include <vector>

double randomval(double lower, double upper)
{
	return lower + (double)rand()/RAND_MAX * (upper - lower);
}

void randomFill(std::vector<double> & out, double lower, double upper)
{
	// std::vector<double> out(dim);
	for (int i = 0; i < out.size(); ++i)
	{
		out[i] = randomval(lower,upper);
	}
	// return out;
}


double ale(double a, double b)
{
	if (a > b)
	{
		return a + log(1.0 + exp(b - a));
	}
	else
	{
		return b + log(1.0 + exp(a - b));
	}
}


double fillFromDeltas(std::vector<double> &x, std::vector<int> & binIds, const std::vector<std::vector<double>> & accumulated_deltas)
{
	double weight = 0;
	int bins = accumulated_deltas[0].size()-1;
	for (int i = 0; i < x.size(); ++i)
	{
		int b = rand() % bins;
		binIds[i] = b;

		double top = accumulated_deltas[i][b+1];
		double bottom = accumulated_deltas[i][b];
		x[i] = randomval(bottom,top);

		weight += log(top-bottom);
	}
	return weight + x.size()*log(bins);
}
void accumulate(const std::vector<std::vector<double>> & deltas, std::vector<std::vector<double>> & target,double lower, double upper)
{
	for (int i =0; i < deltas.size(); ++i)
	{
		target[i][0] = lower;
		for (int j = 0; j < deltas[i].size(); ++j)
		{
			target[i][j+1] = target[i][j] +deltas[i][j];
		}
		target[i][deltas[i].size()] = upper;
	}
}

double fillFromDists(std::vector<double> &x, std::vector<int> & binIds, const std::vector<std::vector<double>> & hists,double lower, double upper)
{
	double log_weight = 0;
	double delta = (upper - lower)/hists[0].size();
	double logdelta = log(delta);
	// std::cout << "\tdelta = " << delta << std::endl;
	// std::cout << "New fill routine " << std::endl;
	for (int i = 0; i < x.size(); ++i)
	{
		
		double r = randomval(0,1);
		// std::cout << "\tDim " << i << " generated " << r << std::endl;
		double run = 0;
		for (int q = 0; q < hists[i].size(); ++q)
		{
			run += hists[i][q];
			if (run >= r)
			{
				
				double binStart = lower + q*delta;
				binIds[i] = q;
				x[i] = randomval(binStart,binStart + delta);
				log_weight += log(hists[i][q]) - logdelta;
				// std::cout << "\t\tDim " << i << " used q = " << q << " = " << binStart << " to generate " << x[i] << std::endl;
				break;
			}
		}
		// std::cout << r << " generated " << histbinIds[i] << " was used to assign! " << std::endl;
	}
	return log_weight;
}