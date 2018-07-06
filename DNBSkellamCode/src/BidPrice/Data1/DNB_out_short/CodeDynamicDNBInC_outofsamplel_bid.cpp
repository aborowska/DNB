#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

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
		// int T = atoi(argv[1]);
		int Model = atoi(argv[1]); 
		// int M = atoi(argv[2]);
		int Burn = atoi(argv[2]);
		int iThin = atoi(argv[3]);
		int Data = atoi(argv[4]);
		
		string sFile, sFile_out, sPrefix, sType;
		string sPath;
		sPath = "/home/aba228/Documents/DNB_Data/";
		
 		switch (Data){
		case 1:
			sFile = "Data_2008103_9_IBM_Bid.csv";
			sFile_out = "Data_20081010_10_IBM_Bid.csv";
			sPrefix = "2008_IBM_";
		    cout << "Data: Bid IMB 2008" << endl;
			break;
		case 2:
			sFile = "Data_2008103_9_KO_Bid.csv";
			sFile_out = "Data_20081010_10_KO_Bid.csv";
			sPrefix = "2008_KO_";
		    cout << "Data: Bid KO 2008" << endl;
			break;
		case 3:
			sFile = "Data_2008103_9_JPM_Bid.csv";
			sFile_out = "Data_20081010_10_JPM_Bid.csv";	
			sPrefix = "2008_JPM_";			
		    cout << "Data: Bid JPM 2008" << endl;
			break;			
		case 4:		
			sFile = "Data_20100423_29_IBM_Bid.csv";
			sFile_out = "Data_20100430_30_IBM_Bid.csv";
			sPrefix = "2010_IBM_";			
		    cout << "Data: Bid IMB 2010" << endl;
			break;			
		case 5:
			sFile = "Data_20100423_29_KO_Bid.csv";
			sFile_out = "Data_20100430_30_KO_Bid.csv";
			sPrefix = "2010_KO_";			
			cout << "Data: Bid KO 2010" << endl;
			break;		
		case 6:
			sFile = "Data_20100423_29_JPM_Bid.csv";
			sFile_out = "Data_20100430_30_JPM_Bid.csv";	
			sPrefix = "2010_JPM_";			
		    cout << "Data: Bid JPM 2010" << endl;
			break;			
		}
		
		
		
		
		if (Model == 1)
		{
			printf("*** Skellam Model Out-of-sample Bid ***\n");
			// printf("M = %i \n", M);
			// printf("BurnIn = %i \n", Burn);
					 
			sType = "Sk";
			/* Parameters*/ 
		//
		//	/* Initialize the MCMC class*/
		//	MCMC NormalMCMC(vNormalParam);
//			MCMC NormalMCMC;
		////	NormalMCMC.SimulateSkellamData(20000);
		////	NormalMCMC.ExportSimulatedData();
		//	NormalMCMC.ImportData("Data_20100423_29_IBM.csv");
		////	NormalMCMC.ImportData("Data_2008103_9_XRX.csv");
		////	NormalMCMC.ImportData("InNormalData.csv");
//			NormalMCMC.ImportData2(sPath, sFile);	
		//
		//	NormalMCMC.EstimateNormal(100000,20000);
//			NormalMCMC.EstimateNormal(M,Burn);
		
//			NormalMCMC.ExportEstimationResults();
		//
		////	cout<<"Estimation Done!"<<endl;
		////	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Skellam","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
		//	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Skellam","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");
		////	OutOfSample Skellam("EstimationResults.csv",20000, 80,"Skellam","InNormalData.csv", "OutNormalData.csv");
			OutOfSample Skellam("EstimationResults.csv",Burn, iThin,sType,sFile,sFile_out,sPath);
		//OutOfSample::OutOfSample(string sFileParam,unsigned int iBurnIn, unsigned int iThin, string sType, string sFileDataIn,string sFileDataOut )
			Skellam.ForecastSkellam(5000); //ForecastSkellam(int iNumOfParticles)
			Skellam.GoodnessOfFitSkellam(1000); //GoodnessOfFitSkellam(int iNumOfParticles)
		}
		else
		{
			/* Parameters*/
			// vector<double> vDNBParam;
				
			printf("*** DNB Model Out-of-sample Bid ***\n");
			// printf("M = %i \n", M);
			// printf("BurnIn = %i \n", Burn);
		 
			sType = "DNB";

			/* Parameters*/ 
		//
		//	/* Initialize the MCMC class*/
		//	MCMC tMCMC(vtParam);
//			MCMC tMCMC;
		//
		////	tMCMC.SimulateDNBData(40000);
		////	tMCMC.ExportSimulatedData();
		////  tMCMC.ImportData("Data_20100423_29_IBM.csv");
		//	tMCMC.ImportData("Data_2008103_9_IBM.csv");
//			tMCMC.ImportData2(sPath, sFile);	
			// tMCMC.ImportData("/home/aba228/Documents/DNB_Data/Data_20100423_29_KO.csv");

		//	tMCMC.EstimateT(100000,20000);
//			tMCMC.EstimateT(M,Burn);
//			tMCMC.ExportEstimationResults();
		
		
		
		
		//	OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
			// OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");
			
			
			OutOfSample DNB("EstimationResults.csv",Burn, iThin,sType,sFile,sFile_out,sPath);
			// OutOfSample DNB("EstimationResults.csv",20000, 80,"DNB","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");

			DNB.ForecastDNB(2000);
			DNB.GoodnessOfFitDNB(1000);
		}
	}
}