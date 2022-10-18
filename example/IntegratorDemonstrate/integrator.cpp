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
std::string fileOutput = "test.eps";
void PreparePlotter(std::vector<int> & dims)
{
	// gp.WindowSize(1700,900);
	gp.SetTerminal("eps");
	gp.SetOutput(fileOutput);
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
	Test(std::string n, double (*func)(int,double,double),bool prevCol,JSL::LineType type)
	{
		Name = n;
		Function = func;
		PrevColour = prevCol;
		Style = type;
	}
	void Prepare(const std::vector<int> & resolutions, const std::vector<int> & amounts)
	{
		int n = resolutions.size();
		Results.resize(n);
		for (int i =0; i < n; ++i)
		{
			Results[i].resize(amounts[i]);
		}

	}
	std::string Name;
	std::vector<std::vector<result>> Results;
	virtual void PerformTest(int loops, int res,int resID,double lower, double upper)
	{
		if (Name.find("GAI") != std::string::npos && res > 1e5)
		{
			Results[resID] = Results[resID-1];
			return;
		}
		auto now = std::chrono::system_clock::now();
		double mVal = 0;
		for (int i = 0; i < loops; ++i)
		{
			
			double v = Function(res,lower,upper);
			mVal += v;
			Results[resID][i].result = v;

		}
		std::chrono::duration<double,std::ratio<1,1000>> duration = std::chrono::system_clock::now() - now;
		for (int i = 0; i < loops; ++i)
		{
			Results[resID][i].time = duration.count()/loops;
		}
		
	}
	typedef double (*functor)(int,double,double); 
	functor Function;
	bool PrevColour;
	JSL::LineType Style = JSL::Solid;
};
struct VegasTest : Test
{
	VegasTest(){}
	VegasTest(std::string n, double (*func)(int,int,int,double,double),int depth,int res,bool prevCol,JSL::LineType type) : Test()
	{
		Name = n;
		vFunction = func;
		vegasDepth = depth;
		vegasResolution = res;
		PrevColour = prevCol;
		Style = type;
	}
	void PerformTest(int loops, int res,int resID,double lower, double upper)
	{
		auto now = std::chrono::system_clock::now();
		// std::vector<result> resultList(loops);
		for (int i = 0; i < loops; ++i)
		{
			
			double v = vFunction(res,vegasDepth,vegasResolution,lower,upper);
			
			Results[resID][i].result = v;

		}
		std::chrono::duration<double,std::ratio<1,1000>> duration = std::chrono::system_clock::now() - now;
		for (int i = 0; i < loops; ++i)
		{
			Results[resID][i].time = duration.count()/loops;
		}
		// Results[resID] = resultList;
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
			double v = 0;
			double t= 0;
			int S = tests[i]->Results[j].size();
			std::vector<double> vals(S);
			double power = 2;
			for (int k = 0; k < S; ++k)
			{
				v += pow(abs(tests[i]->Results[j][k].result - trueVal),power);
				t += tests[i]->Results[j][k].time;
			}	
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
			}
		}
		int pS = 3;
		std::string name = tests[i]->Name;


		if (name == "GAI-Full")
		{
			gp.SetAxis(0,dimID);
			gp.Plot(res,errors,lp::PenSize(pS),lp::Legend(tests[i]->Name));
			gp.SetAxis(1,dimID);
			double mT = JSL::Vector(times).Sum()/times.size();
			double mE = JSL::Vector(errors).Sum()/errors.size();
			gp.Scatter(std::vector<double>{mT},std::vector<double>{mE},lp::ScatterType(JSL::FilledDelta),lp::PenSize(1));
		}
		else
		{
			if (tests[i]->PrevColour && i > 0)
			{
				gp.SetAxis(0,dimID);
				gp.Plot(res,errors,lp::PenSize(pS),lp::Legend(tests[i]->Name),lp::PenType(tests[i]->Style),lp::Colour("hold"));
				gp.SetAxis(1,dimID);
				gp.Plot(times,errors,lp::PenSize(pS),lp::Legend(tests[i]->Name),lp::Colour("hold"),lp::PenType(tests[i]->Style));
			}
			else
			{
				gp.SetAxis(0,dimID);
				gp.Plot(res,errors,lp::PenSize(pS),lp::Legend(tests[i]->Name),lp::PenType(tests[i]->Style));
				gp.SetAxis(1,dimID);
				gp.Plot(times,errors,lp::PenSize(pS),lp::Legend(tests[i]->Name),lp::PenType(tests[i]->Style));
			}
		}
	}

	smallestNonBoring = pow(10,floor(log10(smallestNonBoring)));
	largest = pow(10,ceil(log10(largest)));
	if (smallestNonBoring < 1e-9)
	{
		smallestNonBoring = 1e-9;
	}
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
		gp.SetLegendColumns(2);
		gp.SetFontSize(JSL::Fonts::Legend,5);
	}
	gp.SetXLabel("");
	gp.SetXLabel("Integrator Resolution");
	// if (dimID == dims.size()/2)
	// {
	// 	gp.SetXLabel("Integrator Resolution");
	// }
	gp.SetFontSize(JSL::Fonts::Global,7);
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

}


const double lower = -100;
const double upper = 100;
const int minLoops = 5;
const int maxLoops = 30;
void TestBlock(const std::vector<Test *> & tests, const std::vector<int> & resolutions,const std::vector<int> & amounts,int block, int nBlocks)
{
	auto now = std::chrono::system_clock::now();
	int testCount = tests.size();
	const int resDim = resolutions.size();

	for (int j = 0; j < resDim; ++j)
	{
		for (int i = block; i < testCount; i+=nBlocks)
		{
			tests[i]->PerformTest(amounts[j],resolutions[j],j,lower,upper);
		}
	}
	std::chrono::duration<double,std::ratio<1,1>> duration = std::chrono::system_clock::now() - now;
	std::cout << "\tBlock " << block << " completed in " << JSL::FormatDuration(duration.count()) << std::endl;
}

int main(int argc, char** argv)
{
	auto now = std::chrono::system_clock::now();
	std::cout << "Integrator test initialised" << std::endl;
	JSL::Argument<int> dimInput(1,"d",argc,argv);
	// PreparePlotter
	JSL::Argument<int> seed(time(NULL),"s",argc,argv);
	JSL::Argument<int> testBins(10,"b",argc,argv);
	JSL::Argument<int> testDepth(1,"u",argc,argv);
	JSL::Argument<std::string> file("test.eps","o",argc,argv);
	JSL::Argument<int> nThreads(3,"n",argc,argv);
	const int threadCount = nThreads;
	fileOutput = file;
	srand(seed);

	rand();

	std::vector<int> dims = {10,20,30};
	PreparePlotter(dims);
	int resdim =20;
	int start = 1e1;
	int end = 1e5;
	JSL::Vector res = JSL::Vector::logintspace(start,end,resdim);
	
	resdim = res.Size();
	std::vector<int> resolutions = res;
	std::vector<int> amounts(resdim);
	for (int i = 0; i < resdim; ++i)
	{
		int amount = std::max(minLoops,std::min(maxLoops,(int)(resolutions[resdim-1]/resolutions[i])));
		amounts[i] = amount;
	}
	int dim = dimInput;

	int amount;
	JSL::ProgressBar<2> pb(dims.size(),resdim);
	
	std::vector<std::thread> threads(threadCount-1);

	std::cout << "Thread vector initialised to size " << threadCount << std::endl;
	Test MCI("MCI",&test_MCI,false,JSL::Solid);
	Test LMCI("LMCI",*test_LMCI,true,JSL::Dash);
	Test RGI("RGI",test_RGI,false,JSL::Solid);
	Test LRGI("LRGI",test_LRGI,true,JSL::Dash);
	Test GAI("GAI-Function",test_GAI,true,JSL::Dash);
	Test eGAI("GAI-Full",&exact_GAI,false,JSL::Solid);
	Test bGAI("GAI-BASIC",basic_GAI,false,JSL::Solid);
	int mciD = 5;
	int mciR = 15;
	double fac = 2;
	// VegasTest vMCI("MCI-V",&vegas_MCI,mciD,mciR,false,JSL::Solid);
	VegasTest vLMCI("Log-Vegas",&vegas_LMCI,mciD,mciR,false,JSL::Solid);
	VegasTest vdMCI("Log-Vegas-Half",&vegas_LMCI,mciD,mciR/fac,true,JSL::Dash);
	VegasTest vdLMCI("Log-Vegas-Double",&vegas_LMCI,mciD*fac,mciR*fac,true,JSL::Dotted);
	// VegasTest vvdLMCI("Log-Vegas-Massive",&vegas_LMCI,10,50,false,JSL::Dotted);
	std::cout << "Initialised testers" << threadCount << std::endl;
	// std::vector<Test *> tests = {&RGI,&LRGI,&MCI,&LMCI,&vLMCI,&vdMCI,&vdLMCI,&vvdLMCI,&bGAI,&GAI,&eGAI};

	std::vector<Test *> tests = {&LMCI,&vLMCI,&vdMCI,&vdLMCI,&eGAI};

	for (int q = 0; q < dims.size(); ++q)
	{
		std::cout << "Begin test loop at dimension " << dims[q] << std::endl;
		dim = dims[q];
		means.resize(dim);
		errors.resize(dim);
		randomFill(means,-5,5);
		randomFill(errors,0.02,1);
		
		for (auto t : tests)
		{
			t->Prepare(resolutions,amounts);
		}

		
		for (int b = 0; b < threadCount-1; ++b)
		{
			std::cout << "\tLaunching thread " << b << std::endl;
			threads[b] = std::thread(TestBlock,tests,resolutions,amounts,b,threadCount);
		}
		TestBlock(tests,resolutions,amounts,threadCount-1,threadCount);
		std::cout << "\tLocal thread completed" << std::endl;
		// for (int b = 0; b < threadCount-1; ++b)
		// {
		// 	std::cout << "\tAttempting a join on " << b << std::endl;
		// 	threads[b].join();
		// 	// pb.Update(q,b);
		// }
		int joined = 0;
		while (joined < threadCount - 1)
		{
			for (int b = 0; b < threadCount-1; ++b)
			{
				if (threads[b].joinable())
				{
					threads[b].join();
					++joined;
				}
			}
		}
		// for (int i = 0; i < resdim; ++i)
		// {
		// 	for (int j = 0; j < tests.size(); ++j)
		// 	{
		// 		tests[j]->PerformTest(amounts[i],res[i],lower,upper);
		// 	}
		// 	
		// }

		double trueVal = 0;
		for (int t = 0; t < dim; ++t)
		{
			double d1 = atan((upper-means[t])/errors[t])/M_PI;
			double d2= atan((lower - means[t])/errors[t])/M_PI;
			trueVal += log(d1 - d2);
		}


		plotter(q,dims,res,tests,trueVal);
		// gp.Show();
	}

	gp.Show();
	std::chrono::duration<double,std::ratio<1,1>> duration = std::chrono::system_clock::now() - now;
	std::cout << "Evaluation took " << JSL::FormatDuration(duration.count()) << std::endl;

	return 0;
}