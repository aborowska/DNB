/*
 * MCMC.h
 *
 *  Created on: Mar 17, 2014
 *      Author: istvan
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
	MCMC( vector<double> vParam);
	virtual ~MCMC();

	void SimulateOrderedNormalData(int iNumOfSim);
	void SimulateOrderedTData(int iNumOfSim);
	void ExportSimulatedData();
	void ImportData(string sFile);
	void EstimateNormal( int iNumOfIter, int iBurnIn);
	void EstimateT( int iNumOfIter ,int iBurnIn);
	void ExportEstimationResults();


private:
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
	double* vLambdaTrue;
	double* vSeasonTrue;
	double* vYStarTrue;
	vector<vector<double> > mChain;
	string sInput;
	string sType;
	string sPrefix;

};

#endif /* MCMC_H_ */
