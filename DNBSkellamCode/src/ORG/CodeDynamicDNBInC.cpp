#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
#include "MCMC.h"
#include "OutOfSample.h"



int
main (void)
{





	/* Parameters*/
	// vector<double> vSkellamParam;

	// vSkellamParam.push_back(-1.7);
	// vSkellamParam.push_back(0.97);
	// vSkellamParam.push_back(0.02);
	// vSkellamParam.push_back(0.1);

	// /* Initialize the MCMC class*/
	// MCMC SkellamMCMC(vSkellamParam);

	// SkellamMCMC.SimulateSkellamData(20000);
	// SkellamMCMC.ExportSimulatedData();
	// SkellamMCMC.ImportData("Data_2008103_9_IBM.csv");
	// SkellamMCMC.ImportData("Data_20100423_29_IBM.csv");
	// SkellamMCMC.ImportData("InNormalData.csv");

	// SkellamMCMC.EstimateSkellam(20000,10000);

	// SkellamMCMC.ExportEstaimationResults();
//
//	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Sk","Data_2008103_9_JPM.csv", "Data_20081010_10_JPM.csv");
////	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Sk","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
////	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Sk","InNormalData.csv", "OutNormalData.csv");
//
//
//
//	Skellam.ForecastSkellam(5000);
//	Skellam.GoodnessOfFitSkellam(1000);



	/* Parameters*/
	vector<double> vDNBParam;
	vDNBParam.push_back(-1.7);
	vDNBParam.push_back(0.97);
	vDNBParam.push_back(0.02);
	vDNBParam.push_back(0.2);
	vDNBParam.push_back(10);

	/* Initialize the MCMC class*/
	MCMC DNBMCMC(vDNBParam);
	DNBMCMC.SimulateDNBData(20000);
	DNBMCMC.ExportSimulatedData();
////	DNBMCMC.ImportData("Data_2008103_9_IBM.csv");
////	DNBMCMC.ImportData("Data_20100423_29_KO.csv");
////	DNBMCMC.ImportData("InNormalData.csv");
//
//
//
	DNBMCMC.EstimateDNB(20000,10000);
//
	DNBMCMC.ExportEstimationResults();
//
//
////	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
////	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_2008103_9_IBM.csv", "Data_20081010_10_IBM.csv");
////	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","InNormalData.csv", "OutNormalData.csv");
//
////	DNB.ForecastDNB(5000);
////	DNB.GoodnessOfFitDNB(1000);


}



