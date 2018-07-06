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
	 if ((argc < 8) || (argc >10))
	{
		cout<<"Incorrect number of parameters!"<<endl;
		return -1;
	}
	else
	{
		// int model_Sk = atoi(argv[1]);

		int T = atoi(argv[1]);
		int M = atoi(argv[2]);
		int Burn = atoi(argv[3]);

		// if (model_Sk==1)
		if (argc == 8)
		{
			vector<double> vSkellamParam;	

			printf("*** Skellman Model ***\n");
			printf("T = %i \n", T);
			printf("M = %i \n", M);
			printf("BurnIn = %i \n", Burn);
					
			vSkellamParam.push_back(atof(argv[4]));
			vSkellamParam.push_back(atof(argv[5]));
			vSkellamParam.push_back(atof(argv[6]));
			vSkellamParam.push_back(atof(argv[7]));		
			for (int i=0; i<int(vSkellamParam.size()); i++)
			{
				printf("vSkellamParam[%i] = %6.4f \n", i, vSkellamParam[i]);	
			}
			
			int iSeedNo = 1;	
			/* Initialize the MCMC class*/
			MCMC SkellamMCMC(vSkellamParam, iSeedNo);

			// SkellamMCMC.SimulateSkellamData(50000);
			SkellamMCMC.SimulateSkellamData(T);

			// SkellamMCMC.ExportSimulatedData();
			SkellamMCMC.ExportSimulatedData2("Sk");
			
		//	SkellamMCMC.ImportData("Data_2008103_9_IBM.csv");
		//	SkellamMCMC.ImportData("Data_20100423_29_IBM.csv");
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
			vector<double> vDNBParam;
				
			printf("*** DNB Model ***\n");
			printf("T = %i \n", T);
			printf("M = %i \n", M);
			printf("BurnIn = %i \n", Burn);
		
			vDNBParam.push_back(atof(argv[4]));
			vDNBParam.push_back(atof(argv[5]));
			vDNBParam.push_back(atof(argv[6]));
			vDNBParam.push_back(atof(argv[7]));
			vDNBParam.push_back(atof(argv[8]));
			
			int iSeedNo = atoi(argv[9]);
 			// stringstream ssss;
			// ssss << iSeedNo;
			// cout << ssss.str() << endl;
			// printf("iSeedNo sting %s \n",ssss.str().c_str() ); 
			
			// exit(0);
			for (int i=0; i<int(vDNBParam.size()); i++)
			{
				printf("vDNBParam[%i] = %6.4f \n", i, vDNBParam[i]);	
			}
				
			//	vDNBParam.push_back(-1.7);
			//	vDNBParam.push_back(0.97);
			//	vDNBParam.push_back(0.02);
			//	vDNBParam.push_back(0.2);
			//	vDNBParam.push_back(10);
			
			//	/* Initialize the MCMC class*/
				MCMC DNBMCMC(vDNBParam, iSeedNo);

			//	DNBMCMC.SimulateDNBData(50000);
				DNBMCMC.SimulateDNBData(T);

				// DNBMCMC.ExportSimulatedData();
				DNBMCMC.ExportSimulatedData2("DNB");

			////	DNBMCMC.ImportData("Data_2008103_9_IBM.csv");
			////	DNBMCMC.ImportData("Data_20100423_29_KO.csv");
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