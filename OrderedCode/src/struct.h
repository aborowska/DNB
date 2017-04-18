/*
 * struct.h
 *
 *  Created on: Mar 21, 2014
 *      Author: istvan
 */

#ifndef STRUCT_H_
#define STRUCT_H_
#include <gsl/gsl_linalg.h>

struct AllParam{
	/* Parameters */
	double* dMu;
	double* dPhi;
	double* dSigma2;
	double* dGamma;
	double* dNu;

	/* Auxiliary parameters*/
	int* vN;
	double* vTau1;
	double* vTau2;
	double* mAuxY;
	double* mAuxH;
	double* vZ1;
	double* vZ2;
	double * vYStar;
	double * vLambda;


	/* State */
	double* vX;
	double* vBeta;
	double* vS;

	/* Priors */
	double* dPriorGammaA;
	double* dPriorGammaB;
	double* dPriorNuA;
	double* dPriorNuB;
	double* dPriorMuMean;
	double* dPriorMuSigma2;
	double* dPriorPhiA;
	double* dPriorPhiB;
	double* dPriorSigmaA;
	double* dPriorSigmaB;
	double* dPriorSeasonalMean;
	double* dPriorSeasonalVar;


	/* Num Of Obs */
	int* iNumOfObs;
	int * iNumOfKnots;
	gsl_matrix* W;

};

struct ObjectiveFunParamGamma{
	double* vNormal;
	double* vIndicator;
	int* iNumOfObs;
	double* dPriorA;
	double* dPriorB;
};

struct ObjectiveFunParamPhiSigma{
	int* vN;
	double* mAuxY;
	double* mAuxH;

	int* iNumOfObs;
	double* dPriorPhiA;
	double* dPriorPhiB;
	double* dPriorSigmaA;
	double* dPriorSigmaB;
};



#endif /* STRUCT_H_ */
