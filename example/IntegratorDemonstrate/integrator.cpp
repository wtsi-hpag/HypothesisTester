#include "../HypothesisTester.h"
#include <iostream>
#include <vector>
#include "JSL.h"
#include "targetFunction.h"
#include <algorithm>
#include "mci.h"
#include "funcs.h"
#include "lmci.h"
#include "gai.h"
// JSL::ProgressBar<2>;



void plotter(const std::vector<double> res,const std::vector<std::vector<double>> vals,double trueVal,std::vector<std::string> names)
{


	JSL::gnuplot gp;
	gp.WindowSize(600,800);
	namespace lp = JSL::LineProperties;
	// gp.Plot(res,mci,lp::Legend("MCI"),lp::PenSize(2));
	// gp.Plot(res,lmci,lp::Legend("LMCI"),lp::PenSize(2));
	// gp.Plot(res,lmci_vegas,lp::Legend("LMCI-VEGAS"),lp::PenSize(2));
	// gp.Plot(res,gai,lp::Legend("GAI_{Exact}"),lp::PenSize(2));
	// gp.Plot(res,gai_optim,lp::Legend("GAI"),lp::PenSize(2));
	// gp.SetLegend(true);
	// gp.SetXLog(true);
	// gp.SetYLog(true);
	// double minVal = *std::min_element(lmci.begin(), lmci.end());

	// gp.SetYRange(std::max(minVal,-400.0),5);
	gp.SetMultiplot(2,1);
	std::vector<double> errors(res.size());

	double smallestNonBoring  = 1e300;
	double largest = -1e300;
	for (int i = 0; i < names.size(); ++i)
	{
		gp.SetAxis(0);
		gp.Plot(res,vals[i],lp::Legend(names[i]),lp::PenSize(2));

		for (int j = 0; j < res.size(); ++j)
		{
			errors[j] = abs((vals[i][j] - trueVal));
		}
		double minVal = *std::min_element(errors.begin(),errors.end());
		double maxVal = *std::max_element(errors.begin(),errors.end());
		if (minVal > 1e-12 && minVal < smallestNonBoring && !std::isinf(minVal))
		{
			smallestNonBoring = minVal;
		}
		if (maxVal > largest && !std::isinf(maxVal))
		{
			largest = maxVal;
		}
		std::cout << minVal << "  " << maxVal << "  " << smallestNonBoring << "  " << largest << std::endl;
		gp.SetAxis(1);
		gp.Plot(res,errors,lp::PenSize(2));
	}

	gp.SetAxis(0);
	gp.SetLegend(true);
	gp.SetGrid(true);
	gp.SetXLog(true);
	// gp.SetGlobalColourMap(vals.size(),{1.0,1.0,0.3},{0,0,1});
	gp.SetXRange(*std::min_element(res.begin(),res.end()),*std::max_element(res.begin(),res.end()));
	gp.SetXLabel("Integrator Resolution");
	gp.SetYLabel("Integrator Value");
	// gp.SetYLog(true);
	gp.SetAxis(1);
	gp.SetGrid(true);
	gp.SetXLog(true);
	gp.SetYLog(true);
	gp.SetYRange(smallestNonBoring,largest);
	gp.SetXRange(*std::min_element(res.begin(),res.end()),*std::max_element(res.begin(),res.end()));
	gp.SetYLabel("Log-Fractional Error");
	gp.SetXLabel("Integrator Resolution");
	gp.Show();

}


int main(int argc, char** argv)
{
	JSL::Argument<int> dim(1,"d",argc,argv);
	JSL::Argument<int> seed(time(NULL),"s",argc,argv);
	JSL::Argument<int> testBins(10,"b",argc,argv);
	JSL::Argument<int> testDepth(1,"u",argc,argv);
	srand(seed);

	rand();
	means.resize(dim);
	errors.resize(dim);
	randomFill(means,-8,8);
	randomFill(errors,0.01,2);

	// std::cout << "mus = " << JSL::Vector(means) << std::endl;
	// std::cout << "sigmas = " << JSL::Vector(errors) << std::endl;

	// int b = testBins;
	// int testR = 1e5;
	// int dth = testDepth;
	// double vm = vegas_MCI(testR,dth,b,-10,10);
	// std::cout << "Vegas Estimate = " << vm << std::endl;
	
	// double lm = test_LMCI(testR,-10,10);
	// std::cout << "LMCI Estimate = " << lm << std::endl;
	// exit(5);

	int resdim =40;
	int start = 1e2;
	int end = 1e7;
	JSL::Vector res = JSL::Vector::logintspace(start,end,resdim);

	std::vector<double> mci(resdim);
	std::vector<double> mci_vegas(resdim);
	std::vector<double> lmci(resdim);
	std::vector<double> gai(resdim,cheat_GAI(means));
	std::vector<double> lmci_vegas(resdim);
	std::vector<double> lmci_vegas_1(resdim);
	std::vector<double> lmci_vegas_2(resdim);
	std::vector<double> lmci_vegas_3(resdim);
	std::vector<double> gai_optim(resdim);
	double lower = -30;
	double upper = 30;
	JSL::ProgressBar pb(resdim);
	for (int i = 0; i < resdim; ++i)
	{
		mci[i] = test_MCI(res[i],lower,upper);
		mci_vegas[i] = vegas_MCI(res[i],5,50,lower,upper);
		lmci[i] =(test_LMCI(res[i],lower,upper));
		gai_optim[i] = (test_GAI(res[i],lower,upper));
		lmci_vegas[i] = vegas_LMCI(res[i],5,50,lower,upper);
		lmci_vegas_1[i] = vegas_LMCI(res[i],10,100,lower,upper);
		// lmci_vegas_2[i] = vegas_LMCI(res[i],10,50,lower,upper);
		// lmci_vegas_3[i] = vegas_LMCI(res[i],10,100,lower,upper);
		pb.Update(i);
	}

	
	plotter(res,{mci,mci_vegas,lmci,lmci_vegas,lmci_vegas_1,gai_optim,gai},log(peak),{"MCI","MCI-VEGAS","LMCI","LMCI-VEGAS 5-50","LMCI-VEGAS 5-100","GAI","GAI_{Exact}"});
}