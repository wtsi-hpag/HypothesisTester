#define GNUPLOT_NO_TIDY
#include "../../HypothesisTester.h"
#include <iostream>
#include <vector>
#include "JSL.h"
#include "targetFunction.h"
#include <algorithm>
#include "mci.h"
#include "funcs.h"
#include "lmci.h"
#include "gai.h"
#include "rgi.h"
// JSL::ProgressBar<2>;


JSL::gnuplot gp;

namespace lp = JSL::LineProperties;

void PreparePlotter(std::vector<int> & dims)
{
	gp.WindowSize(1700,900);
	gp.SetTerminal("eps");
	gp.SetOutput("test.eps");
	gp.SetMultiplot(2,dims.size());
}


struct result
{
	double time;
	double result;
};

struct Test
{
	Test(){};
	Test(std::string n, double (*func)(int,double,double))
	{
		Name = n;
		Function = func;
	}
	std::string Name;
	std::vector<std::vector<result>> Results;
	virtual void PerformTest(int loops, int res,double lower, double upper)
	{
		auto now = std::chrono::system_clock::now();
		std::vector<result> resultList(loops);
		for (int i = 0; i < loops; ++i)
		{
			
			double v = Function(res,lower,upper);
			
			resultList[i].result = v;

		}
		std::chrono::duration<double,std::ratio<1,1000>> duration = std::chrono::system_clock::now() - now;
		for (int i = 0; i < loops; ++i)
		{
			resultList[i].time = duration.count()/loops;
		}
		Results.push_back(resultList);
		
	}
	typedef double (*functor)(int,double,double); 
	functor Function;
};
struct VegasTest : Test
{
	VegasTest(){}
	VegasTest(std::string n, double (*func)(int,int,int,double,double),int depth,int res) : Test()
	{
		Name = n;
		vFunction = func;
		vegasDepth = depth;
		vegasResolution = res;
	}
	void PerformTest(int loops, int res,double lower, double upper)
	{
		auto now = std::chrono::system_clock::now();
		std::vector<result> resultList(loops);
		for (int i = 0; i < loops; ++i)
		{
			
			double v = vFunction(res,vegasDepth,vegasResolution,lower,upper);
			
			resultList[i].result = v;

		}
		std::chrono::duration<double,std::ratio<1,1000>> duration = std::chrono::system_clock::now() - now;
		for (int i = 0; i < loops; ++i)
		{
			resultList[i].time = duration.count()/loops;
		}
		Results.push_back(resultList);
	}
	double vegasDepth;
	double vegasResolution;
	typedef double (*vfunctor)(int,int,int,double,double); 
	vfunctor vFunction;
};

void plotter(int dimID, std::vector<int> & dims, const std::vector<double> res, std::vector<Test *> & tests,double trueVal)
{
	int N = res.size();
	std::vector<double> errors(N,0.0);
	std::vector<double> times(N,0);

	double smallestNonBoring  = 1e300;
	double largest = -1e300;
	double shortestTime = 10;
	double longestTime = -100;
	double cut = 1e-14;
	for (int i = 0; i < tests.size(); ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			// std::cout << i << "  " << j << "  " << JSL::Vector(tests[i]->Results[j][k])
			double v = 0;
			double t= 0;
			int S = tests[i]->Results[j].size();
			std::vector<double> vals(S);
			double power = 2;
			for (int k = 0; k < S; ++k)
			{
				v += pow(abs(tests[i]->Results[j][k].result - trueVal),power);
				t += tests[i]->Results[j][k].time;
				// std::cout << "At res " << res[j] << " el " << k << " had val " << tests[i]->Results[j][k].result << std::endl;
			}	
			// std::cout << "Giving mean " << v/S << std::endl;
			errors[j] = pow(v/S,1.0/power)/dims[dimID] + 1e-15;
			times[j] = t/S;
		}

		double minVal = *std::min_element(errors.begin(),errors.end());
		double maxVal = *std::max_element(errors.begin(),errors.end());
		double minTim = *std::min_element(times.begin(),times.end());
		double maxTim = *std::max_element(times.begin(),times.end());
		if (minVal > cut && !std::isinf(minVal))
		{
			if (minVal < smallestNonBoring)
			{
				smallestNonBoring = minVal;
			}
			if (minTim < shortestTime && maxVal < 1e6)
			{
				shortestTime = minTim;
				// std::cout << "New short with " << tests[i]->Name << std::endl;
			}
		}
		if (!std::isinf(maxVal))
		{
			if (maxVal > largest)
			{
				largest = maxVal;
			}
			if (maxTim > longestTime && minVal > cut && maxVal > cut)
			{
				longestTime = maxTim;
				// std::cout << "New long with " << tests[i]->Name << "  " << minVal << "  " << maxVal << std::endl;
			}
		}
		
		std::string name = tests[i]->Name;
		if (name.find("L") == std::string::npos)
		{
			gp.SetAxis(0,dimID);
			gp.Plot(res,errors,lp::PenSize(2),lp::Legend(tests[i]->Name));

			gp.SetAxis(1,dimID);
			if (name.find("GAI-Full")==std::string::npos)
			{
				gp.Plot(times,errors,lp::PenSize(2));
			}
			else
			{
				gp.Scatter(times,errors,lp::ScatterType(JSL::FilledCircle),lp::PenSize(10));
			}
		}
		else
		{
			gp.SetAxis(0,dimID);
			gp.Plot(res,errors,lp::PenSize(2),lp::Legend(tests[i]->Name),lp::Colour("hold"),lp::PenType(JSL::Dash));

			gp.SetAxis(1,dimID);
			
			gp.Plot(times,errors,lp::PenSize(2),lp::Colour("hold"),lp::PenType(JSL::Dash));
		}
		// std::cout << tests[i]->Name << " has " << JSL::Vector(errors) << std::endl;
	}

	smallestNonBoring = pow(10,floor(log10(smallestNonBoring)));
	largest = pow(10,ceil(log10(largest)));
	if (smallestNonBoring < 1e-9)
	{
		smallestNonBoring = 1e-9;
	}
	// std::cout << "Time = " << shortestTime << "  " << longestTime << "\n\n" << std::endl;
	gp.SetAxis(0,dimID);
	gp.SetTitle("Dimension = " + std::to_string(dims[dimID]));
	gp.SetGrid(true);
	gp.SetXLog(true);
	gp.SetYLog(true);
	gp.SetYRange(smallestNonBoring,largest);
	gp.SetYTicPowerFormat(true);
	gp.SetXTicPowerFormat(true);
	gp.SetYLabel("");
	if (dimID == 0)
	{
		gp.SetYLabel("Per-Dimensional Logarithmic Error");
	}
	if (dimID == dims.size() - 1)
	{
		gp.SetLegend(true);
		gp.SetLegendLocation("bottom right");
	}
	gp.SetXLabel("");
	gp.SetXLabel("Integrator Resolution");
	// if (dimID == dims.size()/2)
	// {
	// 	gp.SetXLabel("Integrator Resolution");
	// }
	gp.SetAxis(1,dimID);
	// gp.SetTitle("Dimension = " + std::to_string(dims[dimID]));
	gp.SetGrid(true);
	gp.SetXLog(true);
	gp.SetYLog(true);
	gp.SetYTicPowerFormat(true);
	gp.SetXTicPowerFormat(true);
	gp.SetYRange(smallestNonBoring,largest);
	gp.SetXRange(shortestTime,longestTime);
	gp.SetYLabel("");
	if (dimID == 0)
	{
		gp.SetYLabel("Per-Dimensional Logarithmic Error");
	}
	gp.SetXLabel("Execution Time (ms)");
	// if (dimID == dims.size()/2)
	// {
	// 	gp.SetXLabel("Execution Time (ms)");
	// }
	if (dimID == dims.size() - 1)
	{
		gp.Show();
	}
}


int main(int argc, char** argv)
{
	JSL::Argument<int> dimInput(1,"d",argc,argv);
	// PreparePlotter
	JSL::Argument<int> seed(time(NULL),"s",argc,argv);
	JSL::Argument<int> testBins(10,"b",argc,argv);
	JSL::Argument<int> testDepth(1,"u",argc,argv);
	srand(seed);

	rand();

	std::vector<int> dims = {1,4,12};
	PreparePlotter(dims);
	int resdim =4;
	int start = 1e2;
	int end = 1e6;
	JSL::Vector res = JSL::Vector::logintspace(start,end,resdim);
	resdim = res.Size();
	
	int dim = dimInput;

	int amount;
	JSL::ProgressBar<2,true,'#',50> pb(dims.size(),resdim);
	
	for (int q = 0; q < dims.size(); ++q)
	{
		dim = dims[q];
		means.resize(dim);
		errors.resize(dim);
		randomFill(means,-0,0);
		randomFill(errors,1,1);
		Test MCI("MCI",&test_MCI);
		Test LMCI("LMCI",*test_LMCI);
		Test RGI("RGI",test_RGI);
		Test LRGI("LRGI",test_LRGI);
		Test GAI("GAI-Function",test_GAI);
		Test eGAI("GAI-Full",&exact_GAI);
		Test bGAI("GAI-BASIC",basic_GAI);
		VegasTest vMCI("MCI-V",&vegas_MCI,5,50);
		VegasTest vLMCI("LMCI-V",&vegas_LMCI,5,50);
		VegasTest vdMCI("MCI-V-Deep",&vegas_MCI,20,1000);
		VegasTest vdLMCI("LMCI-Deep",&vegas_LMCI,20,1000);
		double lower = -100;
		double upper = 100;
		
		std::vector<Test *> tests = {&RGI,&LRGI,&MCI,&LMCI,&vMCI,&vLMCI,&bGAI,&GAI,&eGAI};
		// std::vector<Test *> tests = {&RGI,&LMCI,&bGAI,&GAI,&eGAI};
		// std::cout << "Beginning" << std::endl;
		
		for (int i = 0; i < resdim; ++i)
		{
			//exact_GAI(means)
			int ceil = std::min(500, (int)((double)end/res[i]));
			amount = std::max(10,ceil);
			for (int j = 0; j < tests.size(); ++j)
			{
				tests[j]->PerformTest(amount,res[i],lower,upper);
			}
			pb.Update(q,i);
		}
		// double trueVal = 0;//eGAI.Results[0][0].result; //

		double trueVal = 0;
		for (int q = 0; q < dim; ++q)
		{
			double d1 = atan((upper-means[q])/errors[q])/M_PI;
			double d2= atan((lower - means[q])/errors[q])/M_PI;
			trueVal += log(d1 - d2);
		}


		plotter(q,dims,res,tests,trueVal);
	}


}