/*
 * OutOfSample.h
 *
 *  Created on: Aug 19, 2014
 *      Author: istvan
 */

#ifndef OUTOFSAMPLE_H_
#define OUTOFSAMPLE_H_

#include <iostream>
#include <string>
using namespace std;



class OutOfSample {
public:
	OutOfSample(string sFileParam,unsigned int iBurnIn, unsigned int iThin,string sType, string sFileDataIn, string sFileDataOut,string sPath);
	virtual ~OutOfSample();

	void ForecastSkellam(int iNumOfParticles);
	void ForecastDNB(int iNumOfParticles);
	void GoodnessOfFitSkellam( int iNumOfParticles);
	void GoodnessOfFitDNB( int iNumOfParticles);

private:

	/* Number of observations */
	int iNumOfObsIn;
	int iNumOfObsOut;

	/* Number of parameters */
	int iNumOfParam;

	/* Number of variables */
	int iNumOfVar;
	/* Data */
	double* vDataIn;
	double* vTimesIn;
	double* vPriceIn;
	double* vLogReturnIn;
	double* vDataOut;
	double* vTimesOut;
	double* vPriceOut;
	double* vLogReturnOut;
	double* mParam;
	double* vMean;
	string sPrefix;
};

#endif /* OUTOFSAMPLE_H_ */
