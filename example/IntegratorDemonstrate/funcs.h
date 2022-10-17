#pragma once
#include <vector>


std::string comp;
bool isPrinty = false;
void Status(const std::string & title, int res)
{
	if (isPrinty)
	{
		comp = title;
		std::cout << "\tStarting  " << title << " at resolution " << res << std::endl;
	}
}
void Complete()
{
	if (isPrinty)
		std::cout << "\t\t" << comp << " complete " << std::endl;
}

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


inline double ale(double a, double b)
{
	// if (a > b)
	// {
	// 	return a + log(1.0 + exp(b - a));
	// }
	// else
	// {
	// 	return b + log(1.0 + exp(a - b));
	// }
	return std::max(a,b) + log(1.0 + exp(-abs(b-a)));
}


double fillFromDeltas(std::vector<double> &x, std::vector<int> & binIds, const std::vector<std::vector<double>> & accumulated_deltas,double lower, double upper)
{
	double weight = 0;
	const int bins = accumulated_deltas[0].size()-1;
	const double bruch = (double)bins/(upper - lower);
	for (int i = 0; i < x.size(); ++i)
	{
		int b = rand() % bins;
		

		double top = accumulated_deltas[i][b+1];
		double bottom = accumulated_deltas[i][b];
		x[i] = randomval(bottom,top);

		int trueBin = (x[i] - lower) * bruch;
		binIds[i] = trueBin;
		weight += log(top-bottom);
	}
	// std::cout << "Generated " << JSL::Vector(x) << " with weight " << weight +x.size()*log(bins) << std::endl;
	return -weight - x.size()*log(bins);
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

void generateAccumulator(std::vector<std::vector<double>> & accumulate, std::vector<std::vector<double>> & hists, double lower, double upper)
{
	int bins = hists[0].size();
	double delta = (upper - lower)/bins;

	//accumulate contains edges of 1/bin probability
	
	for (int d = 0; d < accumulate.size(); ++d)
	{
		accumulate[d][0] = lower;
		accumulate[d][bins] = upper;
		double runningProb = 0;
		
		int histLoc = 0;
		for (int b = 1; b < bins; ++b)
		{
			double target = (double)b/bins;
			while (runningProb + hists[d][histLoc] <= target)
			{
				runningProb += hists[d][histLoc];
				++histLoc;
			}
			double edge2= lower + (double)(histLoc+1)/(bins) * (upper - lower);
			double edge1 = lower + (double)(histLoc)/(bins) * (upper - lower);
			double bruch =  (target - runningProb)/(hists[d][histLoc]);
			double newX = edge1 + bruch*(edge2 - edge1);
			accumulate[d][b] = newX;
		}
	}
}

void smooth(std::vector<double> & vec)
{
	auto original = vec;

	for (int i =0; i < original.size(); ++i)
	{
		double w = 1;
		double s = original[i];
		double ww = 1;
		int dist = 1;//std::max(1,(int)original.size()/500);
		for (int q = 1; q <= dist; ++q)
		{
			
			if (i-q >= 0)
			{
				s = ale(original[i-q],s);
				w +=ww;
			}
			if (i+q< original.size())
			{
				s = ale(original[i+q],s);
				w+=ww;
			}
		}
		vec[i] = s - log(w);
	}
}