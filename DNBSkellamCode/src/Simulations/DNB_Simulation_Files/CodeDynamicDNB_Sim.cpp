/*
 * CodeDynamicDNB_Sim.cpp
 *
 * The main code executing a single simulation experiment
 * for the Skellam or the DNB model.
 * 
 * Requires the following provided libraries:
 * 		MCMC.h
 *  	OutOfSample.h
 * Requires the following unprovided libraries:
 * 		gsl (GNU Scientific Library)
 * (easily available for Linux)
 * 
 * Input arguments and their types: 
 * 		int T - the required length of the generated time series 
 *		int M - the number of MCMC iterations
 * 		int Burn - the length of the burn-in period
 * 		double mu - the value of mu parameter (unconditional mean of h_t)
 * 		double phi - the value of phi parameter (persistence of the AR(1) process x_t)
 * 		double sigma2 - the value of sigma2_eta parameter (variance of the noise of the AR(1) process x_t)
 * 		double gamma - the value of gamma parameter (zero-inflation)
 * 		double nu - the value of nu parameter (degrees of freedom for DNB)
 *		int iSeedNo - seed number for data generation (index of the simulation)
 *
 * First created by István Barra on Mar 17, 2014
 * Modified by Agnieszka Borowska on Jan 23, 2018
 * 
 */
 
 
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
		int T = atoi(argv[1]);
		int M = atoi(argv[2]);
		int Burn = atoi(argv[3]);

		if (argc == 9) /* Simulate and estimate the Skellam model*/
		{
			vector<double> vSkellamParam;	

			printf("*** Skellman Model ***\n");
			printf("T = %i \n", T);
			printf("M = %i \n", M);
			printf("BurnIn = %i \n", Burn);
					
			/* Set the DGP parameters */		
			vSkellamParam.push_back(atof(argv[4]));
			vSkellamParam.push_back(atof(argv[5]));
			vSkellamParam.push_back(atof(argv[6]));
			vSkellamParam.push_back(atof(argv[7]));
			
			for (int i=0; i<int(vSkellamParam.size()); i++)
			{
				printf("vSkellamParam[%i] = %6.4f \n", i, vSkellamParam[i]);	
			}
			
			/* iSeedNo - seed number for the data simulation, indexes different simulations */
			int iSeedNo = atoi(argv[8]); 
			
			/* Initialize the MCMC class */
			MCMC SkellamMCMC(vSkellamParam, "Sk", iSeedNo);

			/* Simulate data from the Skellam model */			
			SkellamMCMC.SimulateSkellamData(T);
			SkellamMCMC.ExportSimulatedData2();
			
			/* Estimate the Skellam model */			
			SkellamMCMC.EstimateSkellam(M,Burn);
			SkellamMCMC.ExportEstimationResults();

			/* Initialize the OutOfSample class*
			/* e.g. DataIn = __Sk1_SimulatedData.csv */
			//	OutOfSample SkellamOut("EstimationResults.csv",20000, 80,"Sk","DataIn.csv", "DataOut.csv");

			/* Out-of-sample predictive LL and BIC*/						
			//	SkellamOut.ForecastDNB(5000); // number of particles
			
			/* In-sample goodness of fit based on BIC */			
			//	SkellamOut.GoodnessOfFitDNB(1000); // number of particles
			
		}
		else /* Simulate and estimate the DNB model*/
		{
			/* Parameters*/
			vector<double> vDNBParam;
				
			printf("*** DNB Model ***\n");
			printf("T = %i \n", T);
			printf("M = %i \n", M);
			printf("BurnIn = %i \n", Burn);
		
			/* Set the DGP parameters */			
			vDNBParam.push_back(atof(argv[4]));
			vDNBParam.push_back(atof(argv[5]));
			vDNBParam.push_back(atof(argv[6]));
			vDNBParam.push_back(atof(argv[7]));
			vDNBParam.push_back(atof(argv[8]));

			/* iSeedNo - seed number for the data simulation, indexes different simulations */			
			int iSeedNo = atoi(argv[9]);
			
			for (int i=0; i<int(vDNBParam.size()); i++)
			{
				printf("vDNBParam[%i] = %6.4f \n", i, vDNBParam[i]);	
			}

			
			/* Initialize the MCMC class*/
				MCMC DNBMCMC(vDNBParam, "DNB", iSeedNo);
				
			/* Simulate data from the DNB model */						
				DNBMCMC.SimulateDNBData(T);
				DNBMCMC.ExportSimulatedData2();
			
			/* Estimate the Skellam model */						
				DNBMCMC.EstimateDNB(M,Burn);
				DNBMCMC.ExportEstimationResults();
			
			
			/* Initialize the OutOfSample class */
			//	OutOfSample DNBOut("EstimationResults.csv",20000, 80,"DNB","DataIn.csv", "DataOut.csv");
			/* e.g. DataIn = __DNB1_SimulatedData.csv */
			
			/* Out-of-sample predictive LL and BIC*/			
			//	DNBOut.ForecastDNB(5000); // number of particles
			
			/* In-sample goodness of fit based on BIC */
			//	DNBOut.GoodnessOfFitDNB(1000); // number of particles 
		}
	}
}