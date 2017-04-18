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
//	vector<double> vNormalParam;
//
//	vNormalParam.push_back(-1.5);
//	vNormalParam.push_back(0.98);
//	vNormalParam.push_back(0.01);
//	vNormalParam.push_back(0);
//
//	/* Initialize the MCMC class*/
//	MCMC NormalMCMC(vNormalParam);
//
////	NormalMCMC.SimulateOrderedNormalData(20000);
////	NormalMCMC.ExportSimulatedData();
//	NormalMCMC.ImportData("Data_20100423_29_IBM.csv");
////	NormalMCMC.ImportData("Data_2008103_9_XRX.csv");
////	NormalMCMC.ImportData("InNormalData.csv");
//
//
//	NormalMCMC.EstimateNormal(100000,20000);
//
//	NormalMCMC.ExportEstimationResults();
//
////	cout<<"Estimation Done!"<<endl;
////	OutOfSample OrderedNormal("EstimationResults.csv",20000, 80,"OrdNormal","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
//	OutOfSample OrderedNormal("EstimationResults.csv",20000, 80,"OrdNormal","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");
////	OutOfSample OrderedNormal("EstimationResults.csv",20000, 80,"OrdNormal","InNormalData.csv", "OutNormalData.csv");
//
//	OrderedNormal.ForecastOrdNormal(5000);
//	OrderedNormal.GoodnessOfFitOrdNormal(1000);


	/* Parameters*/
//	vector<double> vtParam;
//
//	vtParam.push_back(0);
//	vtParam.push_back(0.95);
//	vtParam.push_back(0.01);
//	vtParam.push_back(0);
//	vtParam.push_back(10);
//
//	/* Initialize the MCMC class*/
//	MCMC tMCMC(vtParam);
//
////	tMCMC.SimulateOrderedTData(40000);
////	tMCMC.ExportSimulatedData();
////  tMCMC.ImportData("Data_20100423_29_IBM.csv");
//	tMCMC.ImportData("Data_2008103_9_IBM.csv");
//
//
//	tMCMC.EstimateT(100000,20000);
//
//	tMCMC.ExportEstimationResults();
//
//	OutOfSample OrderedT("EstimationResults.csv",20000, 80,"OrdT","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
	OutOfSample OrderedT("EstimationResults.csv",20000, 80,"OrdT","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");

	OrderedT.ForecastOrdT(5000);
	OrderedT.GoodnessOfFitOrdT(1000);


}



