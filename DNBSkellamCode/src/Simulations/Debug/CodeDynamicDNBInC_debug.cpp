#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>


using namespace std;
#include "MCMC.h"
#include "OutOfSample.h"

int
main (int argc, char* argv[])
{	

	/* Parameters*/
	 if ((argc < 5) || (argc >5))
	{
		cout<<"Incorrect number of parameters!"<<endl;
		return -1;
	}
	else
	{
		// int model_Sk = atoi(argv[1]);

		// int T = atoi(argv[1]);
		int Model = atoi(argv[1]); 
		int M = atoi(argv[2]);
		int Burn = atoi(argv[3]);
		int Data = atoi(argv[4]);
		
		string sFile; 
		
 		switch (Data){
		case 1:
			sFile = "DNB_010_SimulatedData.csv";
		    cout << "Gamma true 0.1" << endl;
			break;
		case 2:
			sFile = "DNB_001_SimulatedData.csv";
		    cout << "Gamma true 0.01" << endl;
			break;
		case 3:
			sFile = "DNB_050_SimulatedData.csv";
		    cout << "Gamma true 0.5" << endl;
			break;			
		}
		
		
		
		
		if (Model == 1)
		{
			// vector<double> vSkellamParam;	

			printf("*** Skellman Model - DNB sim data  ***\n");
			printf("M = %i \n", M);
			printf("BurnIn = %i \n", Burn);
					 
			
			//int iSeedNo = 1;	
			/* Initialize the MCMC class*/
			// MCMC SkellamMCMC(vSkellamParam, iSeedNo);
			MCMC SkellamMCMC;
			// SkellamMCMC.SimulateSkellamData(50000);
			//SkellamMCMC.SimulateSkellamData(T);

			// SkellamMCMC.ExportSimulatedData();
			//SkellamMCMC.ExportSimulatedData2("Sk");
			
		//	SkellamMCMC.ImportData("Data_2008103_9_IBM.csv");
//			SkellamMCMC.ImportData(sFile);	
			SkellamMCMC.ImportData2("", sFile);		
			// printf("Data imported: Data_20100423_29_IBM\n");

			
		//	SkellamMCMC.ImportData("InNormalData.csv");

			// SkellamMCMC.EstimateSkellam(100000,20000);
			SkellamMCMC.EstimateSkellam(M,Burn);
			SkellamMCMC.ExportEstimationResults();
		//
		//	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Sk","Data_2008103_9_JPM.csv", "Data_20081010_10_JPM.csv");
		////	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Sk","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
		////	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Sk","InNormalData.csv", "OutNormalData.csv");
		//
		//
		//	Skellam.ForecastSkellam(5000);
		//	Skellam.GoodnessOfFitSkellam(1000);

		}
		else
		{
			/* Parameters*/
			// vector<double> vDNBParam;
				
			printf("*** DNB Model -  DNB sim data ***\n");
			printf("M = %i \n", M);
			printf("BurnIn = %i \n", Burn);
		 
			// int iSeedNo = atoi(argv[3]);
 
			//	/* Initialize the MCMC class*/
				// MCMC DNBMCMC(vDNBParam, iSeedNo);
			MCMC DNBMCMC;

			//	DNBMCMC.SimulateDNBData(50000);
			//	DNBMCMC.SimulateDNBData(T);

			//  DNBMCMC.ExportSimulatedData();
			//	DNBMCMC.ExportSimulatedData2("DNB");

			//  DNBMCMC.ImportData("~/Documents/DNB_Data/Data_2008103_9_IBM.csv");
			// DNBMCMC.ImportData("/home/aba228/Documents/DNB_Data/Data_20100423_29_KO.csv");
	//		DNBMCMC.ImportData(sFile);
			DNBMCMC.ImportData2("", sFile);	
			
			////	DNBMCMC.ImportData("InNormalData.csv");

	
			//	DNBMCMC.EstimateDNB(100000,20000);
			DNBMCMC.EstimateDNB(M,Burn);
			//
			DNBMCMC.ExportEstimationResults();
			
			
			////	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
			////	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_2008103_9_IBM.csv", "Data_20081010_10_IBM.csv");
			////	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","InNormalData.csv", "OutNormalData.csv");
			//
			////	DNB.ForecastDNB(5000);
			////	DNB.GoodnessOfFitDNB(1000);
		}
	}
}