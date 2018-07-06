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
	 if ((argc < 4) || (argc >4))
	{
		cout<<"Incorrect number of parameters!"<<endl;
		return -1;
	}
	else
	{
		// int T = atoi(argv[1]);
		int Model = atoi(argv[1]); 
		// int M = atoi(argv[2]);
		// int Burn = atoi(argv[3]);
		int iThin = atoi(argv[2]);
		int Data = atoi(argv[3]);
		
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
			printf("*** OrdNormal Model Out-of-sample Bid ***\n");
			// printf("M = %i \n", M);
			// printf("BurnIn = %i \n", Burn);
					 
			sType = "OrdNormal";
			/* Parameters*/ 
		//
		//	/* Initialize the MCMC class*/
		//	MCMC NormalMCMC(vNormalParam);
//			MCMC NormalMCMC;
		////	NormalMCMC.SimulateOrderedNormalData(20000);
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
		////	OutOfSample OrderedNormal("EstimationResults.csv",20000, 80,"OrdNormal","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
		//	OutOfSample OrderedNormal("EstimationResults.csv",20000, 80,"OrdNormal","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");
		////	OutOfSample OrderedNormal("EstimationResults.csv",20000, 80,"OrdNormal","InNormalData.csv", "OutNormalData.csv");
			OutOfSample OrderedNormal("EstimationResults.csv",20000, iThin,sType,sFile,sFile_out,sPath);
		//OutOfSample::OutOfSample(string sFileParam,unsigned int iBurnIn, unsigned int iThin, string sType, string sFileDataIn,string sFileDataOut )
			OrderedNormal.ForecastOrdNormal(5000); //ForecastOrdNormal(int iNumOfParticles)
			OrderedNormal.GoodnessOfFitOrdNormal(1000); //GoodnessOfFitOrdNormal(int iNumOfParticles)
		}
		else
		{
			/* Parameters*/
			// vector<double> vDNBParam;
				
			printf("*** OrdT Model Out-of-sample Bid ***\n");
			// printf("M = %i \n", M);
			// printf("BurnIn = %i \n", Burn);
		 
			sType = "OrdT";

			/* Parameters*/ 
		//
		//	/* Initialize the MCMC class*/
		//	MCMC tMCMC(vtParam);
//			MCMC tMCMC;
		//
		////	tMCMC.SimulateOrderedTData(40000);
		////	tMCMC.ExportSimulatedData();
		////  tMCMC.ImportData("Data_20100423_29_IBM.csv");
		//	tMCMC.ImportData("Data_2008103_9_IBM.csv");
//			tMCMC.ImportData2(sPath, sFile);	
			// tMCMC.ImportData("/home/aba228/Documents/DNB_Data/Data_20100423_29_KO.csv");

		//	tMCMC.EstimateT(100000,20000);
//			tMCMC.EstimateT(M,Burn);
//			tMCMC.ExportEstimationResults();
		
		
		
		
		//	OutOfSample OrderedT("EstimationResults.csv",20000, 80,"OrdT","Data_20100423_29_XRX.csv", "Data_20100430_30_XRX.csv");
			// OutOfSample OrderedT("EstimationResults.csv",20000, 80,"OrdT","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");
			
			
			OutOfSample OrderedT("EstimationResults.csv",20000, iThin,sType,sFile,sFile_out,sPath);
			// OutOfSample OrderedT("EstimationResults.csv",20000, 80,"OrdT","Data_2008103_9_XRX.csv", "Data_20081010_10_XRX.csv");

			OrderedT.ForecastOrdT(5000);
			OrderedT.GoodnessOfFitOrdT(1000);
		}
	}
}