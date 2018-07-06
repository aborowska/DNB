/*
 * MCMC.h
 *	
 * Definition of the MCMC class (objects and methods)
 * 
 * First created by István Barra on Mar 17, 2014
 * Modified by Agnieszka Borowska on Jan 23, 2018
 * 
 */

#ifndef MCMC_H_
#define MCMC_H_

#include <vector>
#include <iostream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;


class MCMC {
public:
	MCMC();
	MCMC( vector<double> vParam,  string sType, int iSeedNo);
	virtual ~MCMC();

	void SimulateSkellamData(int iNumOfSim);
 	void SimulateDNBData(int iNumOfSim);
	void ExportSimulatedData();
	void ExportSimulatedData2();
	void ImportData(string sFile);
	void ImportData2(string sPath, string sFile);	
	void EstimateSkellam( int iNumOfIter, int iBurnIn);
	void EstimateDNB(int iNumOfIter, int iBurnIn);
	void ExportEstimationResults();
	void BootstrapFilterSkellam( vector<double> vParam,int iNumOfParticles);
	void BootstrapFilterDNB( vector<double> vParam,int iNumOfParticles);


private:	
	/* Seed No */
	int iSeedNo;
	
	/* Parameters*/
	double dMuTrue;
	double dPhiTrue;
	double dSigma2True;
	double dGammaTrue;
	double dNuTrue;

	/* Number of observations */
	int iNumOfObs;

	/* Number of Parameters */
	 int iNumOfParam;

	/* Data */
	double* vData;
	double* vTimes;
	double* vXTrue;
	double* vNTrue;
	double* vTau1True;
	double* vTau2True;
	double* vSeasonTrue;
	vector<vector<double> > mChain;
	string sInput;
	string sType;
	string sPrefix;

};

#endif /* MCMC_H_ */
