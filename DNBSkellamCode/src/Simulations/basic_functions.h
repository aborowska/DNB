/*
 * basic_functions.h
 *
 *  Created on: Mar 18, 2014
 *      Author: istvan
 */

#ifndef BASIC_FUNCTIONS_H_
#define BASIC_FUNCTIONS_H_
#include <sys/stat.h>
#include "struct.h"
#include <math.h>

#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <iomanip>


#define PI 3.141592653589793238462643383279502884197169399375105820974944


int isthisnan(double x){return x!=x;}
int isthisinf(double x){return (!isthisnan(x) && isthisnan(x-x)) | (fabs(x)>DBL_MAX );}
int isthisfinite(double x){return !(isthisnan(x) || isthisinf(x));}


double HyperGeo(double a, double b, double c, double z)
{
	double dTol=0.0000000000000000001;
	double S=1;
	double C=1;
	double dCrit1=1;
	double dCrit2=1;
	double dCrit3=1;

	double SN=1;
	double CN=1;
	double CN1=1;
	double SN1=1;
	double SN2=1;


	int i=0;
	while((dCrit1>dTol)| (dCrit2>dTol)|(dCrit3>dTol)|(i<3))
	{
		if(i==0)
		{
			SN2=S;
		}
		else if(i==1)
		{
			SN1=S;
			CN1=C;

		}
		else if(i==2)
		{
			SN=S;
			CN=C;
		}
		else{
			SN2=SN1;
			SN1=SN;
			SN=S;

			CN1=CN;
			CN=C;

		}

		C=C*((a+i)*(b+i)/(c+i))*(z/(i+1));
		dCrit1=fabs(C)/fabs(S);
		dCrit2=fabs(CN)/fabs(SN1);
		dCrit3=fabs(CN1)/fabs(SN2);

//		cout<<"Crit 1 "<<dCrit1<<endl;
//		cout<<"Crit 2 "<<dCrit2<<endl;
//		cout<<"Crit 3 "<<dCrit3<<endl;
		S=S+C;
		i+=1;
		cout<< i<<endl;
	}
	return S;
}

void CreatPrefix(string sInput,string sType, string *sPrefix )
{
	string sSplit[4];
	string token;
	int i=0;
	while(token != sInput){
	  token = sInput.substr(0,sInput.find_first_of("_"));
	  sInput = sInput.substr(sInput.find_first_of("_") + 1);
	  sSplit[i]=token.c_str();
	  i+=1;
	}

	cout<<sSplit[1].substr(0,4)+"_"+sSplit[3].substr(0,sSplit[3].find(".")) + "_"+sType+"_" <<endl;
	sPrefix[0]=sSplit[1].substr(0,4)+"_"+sSplit[3].substr(0,sSplit[3].find("."))+"_"+sType+"_";
	cout<<sPrefix[0]<<endl;
}


void WriteOutDoubleArray(double * mData, int iRows, int iCols, string  sFile )
{
	mkdir("OutPut",S_IRWXU|S_IRGRP|S_IXGRP);
	fstream ofsData;
	ofsData.open(("OutPut/"+sFile).c_str(), ios::out);
	if(!ofsData.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	for(int i=0; i<iRows; i++)
	{

		for(int j=0; j<iCols ; j++)
		{
			if(j==iCols-1)
			{
				ofsData << mData[i*iCols+j]  << endl;
			}
			else
			{
				ofsData << mData[i*iCols+j]  << ",";
			}
		}
	}
	ofsData.close();
}

void WriteOutIntArray(int * mData,  int iRows, int iCols, string  sFile )
{
	mkdir("OutPut",S_IRWXU|S_IRGRP|S_IXGRP);
	fstream ofsData;
	ofsData.open(("OutPut/"+sFile).c_str(), ios::out);
	if(!ofsData.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	for(int i=0; i<iRows; i++)
	{

		for(int j=0; j<iCols ; j++)
		{
			if(j==iCols-1)
			{
				ofsData << mData[i*iCols+j]  << endl;
			}
			else
			{
				ofsData << mData[i*iCols+j]  << ",";
			}
		}
	}
	ofsData.close();
}

/*
 * Zero inflated zero mean Skellam pdf
 *
 * iX		 	integer argument
 * dGamma 		zero inflation parameter should be between 0 and 1
 * dLambda		intensity parameter
 *
 */
double ZeroSkellamPdf(const int iX, const double dGamma, const double dLambda){

	double dPdf;
	double d2Lambda;
	d2Lambda=2*dLambda;
	if(iX==0)
	{
//		dPdf= dGamma + (1-dGamma)*exp(-d2Lambda)*gsl_sf_bessel_In(0, d2Lambda );
		dPdf= dGamma + (1-dGamma)*gsl_sf_bessel_In_scaled(0, d2Lambda );
	}
	else
	{
		unsigned int iAbsX;
		iAbsX=abs(iX);

//		dPdf=(1-dGamma)*exp(-d2Lambda)*gsl_sf_bessel_In(iAbsX,  d2Lambda);
		dPdf=(1-dGamma)*gsl_sf_bessel_In_scaled(iAbsX,  d2Lambda);
	}

	return dPdf;
}

void SimulateSeasonal( int iNumOfSimPerDay, int iNumOfDays, double* vTimes, double * vSeason)
{
	int iTimeInSec=23400;
	double dDuration=(double)iTimeInSec/iNumOfSimPerDay;

	double *vUnnormalizedSeason=new double[iNumOfSimPerDay*iNumOfDays];
	double *dDailyAvgSeason=new double;



	/* Unnormalized sesonal pattern */
	for(int d=0; d<iNumOfDays; d++)
	{
		dDailyAvgSeason[0]=0;
		for(int i=0; i<iNumOfSimPerDay; i++)
		{
			vTimes[i+d*iNumOfSimPerDay]=(i+1)*dDuration;
			vUnnormalizedSeason[i+d*iNumOfSimPerDay]=(vTimes[i]-10800)*(vTimes[i]-10800)/(58320000);
			dDailyAvgSeason[0]=dDailyAvgSeason[0]+vUnnormalizedSeason[i];
		}

		dDailyAvgSeason[0]=dDailyAvgSeason[0]/iNumOfSimPerDay;

		/* Zero mean normalized seasonal pattern */
		for(int i=0; i<iNumOfSimPerDay; i++)
		{
			vSeason[i+d*iNumOfSimPerDay]=vUnnormalizedSeason[i+d*iNumOfSimPerDay]-dDailyAvgSeason[0];
//			vSeason[i+d*iNumOfSimPerDay]=0;
		}
	}


	/* Saving results */
	string sSeasonFile="vSeason.csv";
	WriteOutDoubleArray(vSeason, iNumOfSimPerDay*iNumOfDays, 1, sSeasonFile);

	string sTimeFile="vTime.csv";
	WriteOutDoubleArray(vTimes, iNumOfSimPerDay*iNumOfDays, 1, sTimeFile);

	delete [] vUnnormalizedSeason;
	delete dDailyAvgSeason;
}


void SimulateSeasonal_str(string sType, int iNumOfSimPerDay, int iNumOfDays, double* vTimes, double * vSeason)
{
	string sPrefix;
	string sInput;
	CreatPrefix(sInput,sType, &sPrefix );
		
	int iTimeInSec=23400;
	double dDuration=(double)iTimeInSec/iNumOfSimPerDay;

	double *vUnnormalizedSeason=new double[iNumOfSimPerDay*iNumOfDays];
	double *dDailyAvgSeason=new double;

	/* Unnormalized sesonal pattern */
	for(int d=0; d<iNumOfDays; d++)
	{
		dDailyAvgSeason[0]=0;
		for(int i=0; i<iNumOfSimPerDay; i++)
		{
			vTimes[i+d*iNumOfSimPerDay]=(i+1)*dDuration;
			vUnnormalizedSeason[i+d*iNumOfSimPerDay]=(vTimes[i]-10800)*(vTimes[i]-10800)/(58320000);
			dDailyAvgSeason[0]=dDailyAvgSeason[0]+vUnnormalizedSeason[i];
		}

		dDailyAvgSeason[0]=dDailyAvgSeason[0]/iNumOfSimPerDay;

		/* Zero mean normalized seasonal pattern */
		for(int i=0; i<iNumOfSimPerDay; i++)
		{
			vSeason[i+d*iNumOfSimPerDay]=vUnnormalizedSeason[i+d*iNumOfSimPerDay]-dDailyAvgSeason[0];
//			vSeason[i+d*iNumOfSimPerDay]=0;
		}
	}


	/* Saving results */
	string sSeasonFile=(sPrefix+"vSeason.csv").c_str();
	WriteOutDoubleArray(vSeason, iNumOfSimPerDay*iNumOfDays, 1, sSeasonFile);

	string sTimeFile=(sPrefix+"vTime.csv").c_str();
	WriteOutDoubleArray(vTimes, iNumOfSimPerDay*iNumOfDays, 1, sTimeFile);

	delete [] vUnnormalizedSeason;
	delete dDailyAvgSeason;
}

void SimulateSeasonal2( int iNumOfSim, double* vTimes, double * vSeason)
{
	int iTimeInSec=23400;
	double dDuration=(double)iTimeInSec/iNumOfSim;

	double *vUnnormalizedSeason=new double[iNumOfSim];
	double *dAvgSeason=new double;

	dAvgSeason[0]=0;

	/* Unnormalized sesonal pattern */
	for(int i=0; i<iNumOfSim; i++)
	{
		vTimes[i]=(i+1)*dDuration;
		if(i<iNumOfSim/4)
		{
			vUnnormalizedSeason[i]=0;
		}
		else if ((i<iNumOfSim/2) & (i >= iNumOfSim/4))
		{
			vUnnormalizedSeason[i]=0;
		}
		else if ((i<3*iNumOfSim/4) & (i >= iNumOfSim/2))
		{
			vUnnormalizedSeason[i]=0;
		}
		else
		{
			vUnnormalizedSeason[i]=0;
		}

		dAvgSeason[0]=dAvgSeason[0]+vUnnormalizedSeason[i];
	}

	dAvgSeason[0]=dAvgSeason[0]/iNumOfSim;

	/* Zero mean normalized sesonal pattern */
	for(int i=0; i<iNumOfSim; i++)
	{
		vSeason[i]=vUnnormalizedSeason[i]-dAvgSeason[0];
	}

	/* Saving results */
	string sSeasonFile="vSeason.csv";
	WriteOutDoubleArray(vSeason, iNumOfSim, 1, sSeasonFile);

	string sTimeFile="vTime.csv";
	WriteOutDoubleArray(vTimes, iNumOfSim, 1, sTimeFile);

	delete [] vUnnormalizedSeason;
	delete dAvgSeason;
}

void SimulateSkellam(double dLambda, const gsl_rng * gsl_random_num, int*  iY, int* iN, double* dTau1, double * dTau2)
{
	double dTau=0;
	double dTauPrev=0;

	int iNTemp=0;
	int iYTemp=0;

	int iNPrevTemp=0;
	int iYPrevTemp=0;

	while(dTau < 1)
	{
		iNPrevTemp=iNTemp;
		iYPrevTemp=iYTemp;
		dTauPrev=dTau;
		dTau=dTau+gsl_ran_exponential(gsl_random_num, 1/(2*dLambda));

		iNTemp=iNTemp+1;
		double dU=gsl_rng_uniform(gsl_random_num);
		if(dU<=0.5)
		{
			iYTemp=iYTemp+1;
		}
		else
		{
			iYTemp=iYTemp-1;
		}
	}

	iN[0]=iNPrevTemp;
	iY[0]=iYPrevTemp;
	dTau1[0]=dTau-dTauPrev;
	dTau2[0]=dTauPrev;
	// cout <<"while end"<< endl;
}

double CalcVol(double* vData, int iCurrent ,int iWindow)
{
	double dSS=0;
	double dS=0;
	for(int i=0; i<iWindow; i++)
	{
		dS=dS+vData[iCurrent-i];
		dSS=dSS+vData[iCurrent-i]*vData[iCurrent-i] ;
	}

	return dSS/iWindow-(dS*dS)/(iWindow*iWindow);
}

void CalculateInitialLogVol(double* vData, double * vTime,  int iNumOfObs,  int iWindow,double * vInitialLogVol)
{
	for(int i=iWindow-1; i<iNumOfObs;i++)
	{
		double dVar=CalcVol(vData, i, iWindow);
		if(dVar==0)
		{
			vInitialLogVol[i]=log(0.01);
		}
		else
		{
			vInitialLogVol[i]=log(dVar/2);
		}
	}

	for(int i=0; i<iWindow-1;i++)
	{
		vInitialLogVol[i]=vInitialLogVol[iWindow-1]-(iWindow-1-i)*(vInitialLogVol[2*iWindow-1]-vInitialLogVol[iWindow-1])/iWindow;
	}
}

void CalculateInitialState(struct AllParam sParam, double * vInitialLogVol, int iNumOfObs, int iNumOfKnots,   gsl_matrix * mWTilde  )
{
	gsl_matrix * mYBar=gsl_matrix_alloc (iNumOfObs, 1);
	gsl_matrix * mWBar=gsl_matrix_alloc(iNumOfObs,iNumOfKnots-1);

	double * dY=new double;

	double dMeanLogInt=0;
	for(int i=0; i<iNumOfObs; i++)
	{
		dMeanLogInt=dMeanLogInt+vInitialLogVol[i]/iNumOfObs;

	}

	//int iTempCount=0;
	for(int i=0; i<iNumOfObs; i++)
	{

//			cout << "dMeanLogInt "<<dMeanLogInt << endl;
			dY[0]=vInitialLogVol[i]-dMeanLogInt;


			gsl_matrix_set (mYBar,i,0, dY[0] );
//			cout <<"mYBar "<<gsl_matrix_get (mYBar,i,0)<<endl;
			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,i,j, gsl_matrix_get (mWTilde, i, j));
//				cout <<"mWBar "<<gsl_matrix_get (mWBar,i,j)<<endl;
			}


	}


	/* Variance */
	gsl_matrix * mWW=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVar=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVarInv=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mWBar,0.0, mWW);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			if(j==i)
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j)+sParam.dPriorSeasonalVar[0]);
			}
			else
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j));
			}
		}
	}

	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots-1);
	int s;
	gsl_linalg_LU_decomp (mVarInv, mPermutation, &s);
	gsl_linalg_LU_invert (mVarInv, mPermutation, mVar);



	/* Mean */
	gsl_matrix * mWY=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMean=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMeanTemp=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mYBar,0.0, mWY);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		gsl_matrix_set (mMeanTemp,i,0,gsl_matrix_get(mWY,i,0)+sParam.dPriorSeasonalMean[0]/sParam.dPriorSeasonalVar[0]);
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mVar,mMeanTemp,0.0, mMean);



	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(mMean,j,0)*gsl_matrix_get(mWTilde,i,j);
		}

	}

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vX[i]=vInitialLogVol[i]-dMeanLogInt-sParam.vS[i];
	}

	sParam.dMu[0]=dMeanLogInt;

	gsl_matrix_free (mYBar);
	gsl_matrix_free (mWBar);
	gsl_matrix_free (mWW);
	gsl_matrix_free (mWY);
	gsl_matrix_free (mMean);
	gsl_matrix_free (mMeanTemp);
	gsl_matrix_free (mVar);
	gsl_matrix_free (mVarInv);

	delete dY;

}


double LogModifiedBesselFirstKindLargeX(int iN, double dX)
{
	int iNumOfK=200;
	double dValue=1;
	for(int i=iNumOfK; i>0;i--)
	{
		dValue=1-dValue*(4*pow((double) iN,2)-pow((double) 2*i-1,2))/(8*dX*i);
	}
//	dValue=dValue*exp(dX)/(sqrt(2*PI*dX));
	dValue=dX-0.5*log(2*PI*dX)+log(dValue);
	return dValue;
}

void CaculateSkellam(const double * vData,  int iNumOfObs,  struct AllParam sParam, double * vZeroSkellamDens,
		double * vIndicator,double * vSkellam,  double * vLogBessel  )
{
	double* dLambda=new double;

	for(int i=0; i<iNumOfObs; i++)
	{
		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);
		
//		cout << "Skellam " << i << endl;
//		cout << "iAbsX " << abs(vData[i])<< endl;
//		cout << "-2*dLambda[0] " << -2*dLambda[0]<< endl;
//		cout << "sParam.dSigma[0] " << sParam.dSigma[0]<< endl;
//		cout << "sParam.vX[i] " << sParam.vX[i]<< endl;
		//cout << "not scaled " << exp(-2*dLambda[0])*gsl_sf_bessel_In(abs(vData[i]), 2*dLambda[0])<< endl;
		//cout << "scaled " << gsl_sf_bessel_In_scaled(abs(vData[i]), 2*dLambda[0])<< endl;

		gsl_sf_result result;
		gsl_set_error_handler_off();
		int status = gsl_sf_bessel_In_scaled_e(abs(vData[i]), 2*dLambda[0],&result);
		gsl_set_error_handler (NULL);

		if (status) {
			printf ("error: %s\n", gsl_strerror (status));
			if(abs(vData[i])>70)
			{
				result.val=0;
			}
			else
			{
				result.val=exp(log(-2*dLambda[0])+abs(vData[i]) * log(dLambda[0])-lgamma(abs(vData[i])+1));
			}
		}

		if(1e-50>2*dLambda[0])
		{
			vSkellam[i]=pow(dLambda[0],abs(vData[i]) )/tgamma(abs(vData[i])+1);
		}
		else
		{
			vSkellam[i]=result.val;
		}

		if(vData[i]==0)
		{
			vIndicator[i]=1;
		}
		else
		{
			vIndicator[i]=0;
		}

//		cout << "2*dLambda[0] " <<2*dLambda[0]<<endl;
		if(1e-50>2*dLambda[0])
		{
//			vBessel[i]=pow(dLambda[0],abs(vData[i]) )/tgamma(abs(vData[i])+1);
			vLogBessel[i]=abs(vData[i]) * log(dLambda[0])-lgamma(abs(vData[i])+1);
		}
		else if(700<2*dLambda[0])
		{
			vLogBessel[i]=LogModifiedBesselFirstKindLargeX(abs(vData[i]), 2*dLambda[0]);
		}
		else
		{
			if(result.val==0)
			{
				vLogBessel[i]=-100;
			}
			else{
				vLogBessel[i]=log(result.val);
			}
		}
		vZeroSkellamDens[i]=sParam.dGamma[0]*vIndicator[i]+(1-sParam.dGamma[0])*vSkellam[i] ;
	}
	delete dLambda;
}

/*
 * Zero inflated zero mean negative binomial difference pdf
 *
 * iX		 	integer argument
 * dGamma 		zero inflation parameter should be between 0 and 1
 * dLambda		intensity parameter
 * dNu			"degree of freedom" parameter
 *
 */
double ZeroDNBPdf(const int iX, const double dGamma, const double dLambda, const double dNu){

		double dPdf;
		double dLambdaRatio=dLambda/(dLambda+dNu);
		double dNuRatio=dNu/(dLambda+dNu);

		if(iX==0)
		{
			dPdf= dGamma + (1-dGamma)*pow(dNuRatio,2*dNu)*gsl_sf_hyperg_2F1(dNu,dNu, 1 , pow(dLambdaRatio,2));
		}
		else
		{
			int iAbsX;
			iAbsX=abs(iX);
			dPdf=(1-dGamma)*pow(dNuRatio,2*dNu)*pow(dLambdaRatio,iAbsX)*(gsl_sf_gamma(dNu+iAbsX)/(gsl_sf_gamma(iAbsX+1)*gsl_sf_gamma(dNu)))*gsl_sf_hyperg_2F1(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) );
		}
		return dPdf;
}

int MixtureLogGamma( double n,  int * iNumOfComp,  double * vWeights , double * vMeans ,double *vVariances)
{
	 /*
	  * This routine was initially written by Rudolf Fruehwirth (Institut fuer Hochenergiephysik Oesterreichische Akademie der
	  * Wissenschaften, Nikolsdorfer Gasse 18, 1050 Wien), and later adapted to GMRFLib by H.Rue
	  *
	  * This routine approximate the log-gamma distribution, with density
	  *
	  * exp(-n*x - exp(-x)) / Gamma(n)
	  *
	  * by a mixture of normals
	  *
	  * \sum_{i=1}^N weights[i] * N(x; means[i], variances[i])
	  *
	  * with N components.
	  *
	  */

      const int nrange = 7;
      const int range[7][2] = {
              {1, 4},                                        /* 10 components */
              {5, 19},                                       /* 9 components */
              {20, 49},                                      /* 4 components */
              {50, 439},                                     /* 3 components */
              {440, 1599},                                   /* 2 components */
              {1600, 10000},                                 /* 2 components */
              {10000, 30000}                                 /* 2 components */
      };
      const int ncomp[] = { 10, 9, 4, 3, 2, 2, 2 };
      const double P10[4][10] = {
              {0.00396984425, 0.0396244597, 0.16776747, 0.147036501, 0.125306271,
               0.1014852, 0.103758531, 0.115972617, 0.107065659, 0.0880134482},
              {0.00396860865, 0.039627049, 0.16777003, 0.147037759, 0.125304523,
               0.101481714, 0.103759705, 0.115973128, 0.107065554, 0.0880119305},
              {0.00396904734, 0.0396280642, 0.167770514, 0.14703607, 0.12530043,
               0.10148242, 0.103759287, 0.115974323, 0.107066971, 0.0880128738},
              {0.00396883344, 0.039627504, 0.167771274, 0.147036659, 0.125301189,
               0.101481755, 0.103760036, 0.115974339, 0.107065718, 0.0880126919}
      };
      const double M10[4][10] = {
              {3.55887454, 2.11415904, 0.968124631, 0.51537638, 0.145465449,
               -0.145346445, -0.416660312, -0.689002855, -0.974965634, -1.27310004},
              {2.78807754, 1.84979328, 0.94844169, 0.577673108, 0.223449219,
               -0.0831666379, -0.387174155, -0.69613969, -1.01843553, -1.34112844},
              {2.43454312, 1.7315327, 0.942157786, 0.60208557, 0.251664821,
               -0.0644746918, -0.379817508, -0.696781518, -1.0293035, -1.35705784},
              {2.30484474, 1.71231656, 0.952907078, 0.601128034, 0.252368847,
               -0.059032783, -0.375605704, -0.699071542, -1.03734211, -1.3609072}
      };
      const double V10[4][10] = {
              {2.62603032, 1.21263644, 0.66586521, 0.256650604, 0.120071142,
               0.0649909219, 0.0473513798, 0.046639443, 0.0576541602, 0.0888536903},
              {2.39619753, 1.16995764, 0.688870128, 0.307084756, 0.155644328,
               0.0899360571, 0.0707828448, 0.0751755614, 0.0990773728, 0.15471843},
              {2.16215586, 1.11062998, 0.682294453, 0.324750601, 0.173204837,
               0.108063698, 0.0917073596, 0.100257256, 0.131371692, 0.200024832},
              {1.92939547, 1.00671896, 0.638983371, 0.322852776, 0.18445103,
               0.122217472, 0.106400052, 0.116936918, 0.154113316, 0.233525098}
      };
      const double P9[15][9] = {
              {0.0435820277, 0.167794347, 0.147040722, 0.125310654, 0.10147112,
               0.10376347, 0.115973878, 0.107056197, 0.0880075845},
              {0.0435817033, 0.167795985, 0.1470426, 0.125311016, 0.101470666,
               0.103763084, 0.115972864, 0.107055437, 0.0880066471},
              {0.0435798639, 0.167797087, 0.147042073, 0.125313291, 0.101470979,
               0.103761847, 0.115973234, 0.107054351, 0.0880072753},
              {0.043578895, 0.167797426, 0.147041988, 0.125313875, 0.101470922,
               0.103761581, 0.115973137, 0.107054001, 0.0880081751},
              {0.0435786725, 0.167797743, 0.1470428, 0.125313553, 0.101470946,
               0.103761391, 0.115973188, 0.10705364, 0.0880080663},
              {0.0435779307, 0.167797788, 0.147042734, 0.125314068, 0.101471449,
               0.10376142, 0.115973187, 0.107053473, 0.0880079505},
              {0.043576761, 0.167801375, 0.147042624, 0.125314075, 0.101470546,
               0.103761069, 0.115973226, 0.107051966, 0.0880083593},
              {0.0435771819, 0.167801103, 0.147042441, 0.125313864, 0.101470305,
               0.103761519, 0.11597319, 0.107052417, 0.0880079809},
              {0.0435778469, 0.167800486, 0.147041951, 0.125313914, 0.101470076,
               0.103761707, 0.115973611, 0.107052756, 0.0880076518},
              {0.0435786417, 0.16779926, 0.147042119, 0.125313391, 0.101470554,
               0.103762378, 0.115973792, 0.107052537, 0.0880073289},
              {0.043581505, 0.167797871, 0.147043608, 0.125312521, 0.101469081,
               0.103762173, 0.115973414, 0.107054363, 0.0880054639},
              {0.0435811435, 0.167798952, 0.147043687, 0.125312616, 0.101468918,
               0.103762052, 0.115973417, 0.107053968, 0.0880052462},
              {0.0435812603, 0.167798873, 0.147044518, 0.125312321, 0.101468879,
               0.103761729, 0.115972692, 0.107054049, 0.0880056789},
              {0.0435808733, 0.167799002, 0.147044529, 0.125312675, 0.101468951,
               0.103761472, 0.115972643, 0.107053883, 0.0880059719},
              {0.0435807283, 0.167799231, 0.14704464, 0.12531292, 0.101468814,
               0.103761275, 0.115972628, 0.107053662, 0.088006103}
      };
      const double M9[15][9] = {
              {1.31113348, 0.963928895, 0.659198795, 0.240742429, -0.108844644,
               -0.252087404, -0.6546691, -1.04146524, -1.37874376},
              {1.25919247, 0.957217299, 0.66710982, 0.251658342, -0.125234491,
               -0.240137829, -0.64912733, -1.03921002, -1.37439461},
              {1.21602216, 0.94778507, 0.671484869, 0.265435387, -0.104709908,
               -0.24708343, -0.653441223, -1.04076324, -1.36988994},
              {1.18027937, 0.939725546, 0.67760436, 0.293497817, -0.110879079,
               -0.257696481, -0.655613756, -1.0406543, -1.36465528},
              {1.14996911, 0.934206664, 0.686267712, 0.311595579, -0.112948479,
               -0.274222612, -0.653808807, -1.04092104, -1.35962481},
              {1.12841748, 0.932206841, 0.69102714, 0.319038554, -0.109581301,
               -0.302963892, -0.641448217, -1.03858769, -1.35274157},
              {1.10451126, 0.925180162, 0.689947194, 0.309587296, -0.123979787,
               -0.239246368, -0.658582798, -1.03932069, -1.347407},
              {1.08624068, 0.918081034, 0.697616213, 0.330655882, -0.106424319,
               -0.290644969, -0.644517493, -1.04099153, -1.34370607},
              {1.0671125, 0.915784215, 0.70024231, 0.330800476, -0.125598534,
               -0.244656951, -0.661886313, -1.04447342, -1.33948264},
              {1.05721516, 0.918099637, 0.698999193, 0.325014717, -0.153165358,
               -0.225909041, -0.659788653, -1.03711782, -1.33064663},
              {1.02150943, 0.896206397, 0.702224917, 0.344137939, -0.119895501,
               -0.256590721, -0.641185291, -1.03810889, -1.32943558},
              {1.02508782, 0.902555642, 0.699256309, 0.336391119, -0.121902141,
               -0.242730179, -0.6538063, -1.0385784, -1.32415888},
              {0.997274184, 0.88197491, 0.696155279, 0.3460138, -0.128981232,
               -0.227346713, -0.630322077, -1.03647508, -1.32316505},
              {0.995086849, 0.891409674, 0.70171109, 0.341992158, -0.127906113,
               -0.245952673, -0.638792902, -1.03392281, -1.31486719},
              {0.997741814, 0.892112396, 0.698155553, 0.337731787, -0.122630195,
               -0.240878604, -0.651951415, -1.02898878, -1.3062535}
      };
      const double V9[15][9] = {
              {1.5732832, 0.745075965, 0.340530976, 0.206325108, 0.206977107,
               0.133034557, 0.123981078, 0.155417698, 0.247661591},
              {1.52550277, 0.745216293, 0.347702459, 0.213195645, 0.220928839,
               0.147502243, 0.139478204, 0.17271313, 0.269719569},
              {1.48970429, 0.74910777, 0.35810967, 0.221696291, 0.216470192,
               0.155837875, 0.148481868, 0.185394632, 0.28822907},
              {1.46105103, 0.752441091, 0.365198621, 0.220104509, 0.199190433,
               0.167708126, 0.15761138, 0.197076001, 0.304425302},
              {1.43764551, 0.754190306, 0.367534375, 0.215776065, 0.185257157,
               0.180142183, 0.165402413, 0.206954388, 0.318591695},
              {1.41468216, 0.75198881, 0.368357589, 0.215271168, 0.178178434,
               0.198636491, 0.176790288, 0.218155881, 0.332156859},
              {1.39851898, 0.755429842, 0.377058085, 0.229287048, 0.214645547,
               0.18489307, 0.178139004, 0.226237823, 0.343708183},
              {1.38111403, 0.759024378, 0.379809227, 0.222659694, 0.185443843,
               0.206181273, 0.184773494, 0.231840962, 0.353714302},
              {1.36922993, 0.759197249, 0.381687395, 0.225704876, 0.199623554,
               0.195711194, 0.18270427, 0.236837387, 0.363050264},
              {1.35398708, 0.753650144, 0.381381699, 0.231057971, 0.208319112,
               0.210190241, 0.194957855, 0.249236388, 0.373774124},
              {1.35652837, 0.774748407, 0.400413698, 0.238592235, 0.199467639,
               0.230239828, 0.19924794, 0.251600772, 0.380054821},
              {1.33546695, 0.763749521, 0.396745563, 0.241905327, 0.212176877,
               0.218950701, 0.201882762, 0.257807637, 0.388524892},
              {1.33575722, 0.781739895, 0.417799104, 0.256728889, 0.211230256,
               0.254750255, 0.208700024, 0.26087813, 0.393822372},
              {1.3227643, 0.771070524, 0.406631212, 0.249617029, 0.210958958,
               0.249496089, 0.214362668, 0.270024593, 0.402615529},
              {1.30630549, 0.765952536, 0.407914566, 0.255018833, 0.226289944,
               0.236887588, 0.221124118, 0.280039124, 0.411219814}
      };
      /*
       * n from 20 to 49
       */
#define size_Coeff_p 4
#define size_Coeff_m 5
#define size_Coeff_v 5
      const double Coeff_p3[size_Coeff_p][4] = {
              {-5.644536495326e-009, 7.299190941772e-009, -1.788056445701e-008, 9.259794020097e-009},
              {-1.266992312621e-006, 1.387196986613e-006, -2.391966642312e-006, 1.224613603301e-006},
              {1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000},
              {4.730022618640e+000, 3.672627139064e+000, 4.871566292199e+000, 3.215154075256e+000}
      };
#define size_Coeff_m3 5
      const double Coeff_m3[size_Coeff_m][4] = {
              {4.552797222246e-005, 2.284729919322e-005, -3.900177124794e-005, -2.486737015928e-005},
              {4.638009105861e-002, -1.095058888700e-002, 4.731686443506e-002, 2.978371498898e-002},
              {1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000},
              {9.627160143020e-002, -1.690501643196e-002, -5.610095109269e-001, -4.643825308040e-002},
              {1.143772956136e+000, 1.944583776810e+000, -6.033854619021e+000, -1.105498133467e+000}
      };
#define size_Coeff_v3 5
      const double Coeff_v3[size_Coeff_v][4] = {
              {-2.191015160635e-005, 7.060864706965e-005, 1.823003483481e-004, 1.613752763707e-004},
              {9.939739739229e-002, 1.143203813438e-001, 1.675101633325e-001, 1.943336591437e-001},
              {1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000},
              {9.208564449364e-002, 1.548617268518e-001, 2.735383307281e-001, 2.797653349940e-001},
              {7.148740127686e-001, 2.428636911969e+000, 4.861423133312e+000, 3.840341872065e+000}
      };
      /*
       * n from 50 to 439
       */
      const double Coeff_p4[size_Coeff_p][3] = {
              {-5.639545796991280e-010, 2.651836392450035e-010, -2.384482520627535e-011},
              {4.698743002874532e-007, -1.280380156002802e-007, -1.227680572544847e-007},
              {1.000000000000000e+000, 1.000000000000000e+000, 1.000000000000000e+000},
              {4.730482920811330e+000, 2.093982718501769e+000, 3.214956149674574e+000}
      };
      const double Coeff_m4[size_Coeff_m][3] = {
              {-1.653173201148335e-006, -8.298537364426537e-007, -1.431525987300163e-006},
              {1.036578627632170e-002, 5.017456263052972e-003, 8.386323466104712e-003},
              {1.000000000000000e+000, 1.000000000000000e+000, 1.000000000000000e+000},
              {2.349390607953272e-002, 5.123168011502721e-002, -1.841057020139425e-002},
              {1.432904956685477e+000, 6.453910704667408e+000, -1.410602407670769e+000}
      };
      const double Coeff_v4[size_Coeff_v][3] = {
              {-2.726183914412441e-007, 1.118379212729684e-006, 2.197737873275589e-006},
              {2.788507874891710e-002, 2.433214514397419e-002, 3.186581505796005e-002},
              {1.000000000000000e+000, 1.000000000000000e+000, 1.000000000000000e+000},
              {2.777086294607445e-002, 2.778340896223197e-002, 3.808382220884354e-002},
              {8.369406298984288e-001, 1.489387981224663e+000, 1.958805931276004e+000}
      };
      /*
       * n from 440 to 1599
       */
      const double Coeff_p5[size_Coeff_p][2] = {
              {1.034981541036597e-010, -2.291586556531707e-010},
              {-2.445177000398938e-007, 5.414543692806514e-007},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {1.451229377475864e+000, 3.216167113242079e+000}
      };
      const double Coeff_m5[size_Coeff_m][2] = {
              {-6.578325435644067e-008, -6.292364160498604e-008},
              {1.648723149067166e-003, 1.618047470065775e-003},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {1.594968525045459e-002, -7.091699113800587e-003},
              {5.566082591106806e+000, -2.516741952410371e+000}
      };
      const double Coeff_v5[size_Coeff_v][2] = {
    		  {-2.802162650788337e-009 ,   3.776558110733883e-008},
    		  {4.051503597935380e-003 ,   5.022619018941299e-003},
    		  {1.000000000000000e+000  ,  1.000000000000000e+000},
    		  {4.018981069179972e-003  ,  5.525253413878772e-003},
    		  {9.654061278849895e-001  ,  1.450507513327352e+000},
      };

      /*
       * n from 1600 to 10000
       */
      const double Coeff_p6[size_Coeff_p][2] = {
              {-1.586037487404490e-013, 1.291237745205579e-013},
              {3.575996226727867e-009, -2.911316152726367e-009},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {2.228310599179340e+000, 1.814126328168031e+000}
      };
      const double Coeff_m6[size_Coeff_m][2] = {
              {-2.419956255676409e-009, -2.419411092563945e-009},
              {3.245753451748892e-004, 3.245669014250788e-004},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {1.895335618211674e-003, -2.327930564510444e-003},
              {3.388553853864067e+000, -4.162274939236667e+000}
      };
      const double Coeff_v6[size_Coeff_v][2] = {
              {-6.024563976875348e-011, 5.024317053887777e-010},
              {-6.540694956580495e-004, 8.898044793516080e-004},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {-6.582951415419203e-004, 9.246987493760628e-004},
              {1.006399508694657e+000, 1.149073788967684e+000}
      };
      /*
       * n from 10000 to 30000
       */
      const double Coeff_p7[size_Coeff_p][2] = {
              {-1.663426552872397e-014, 1.354267905471566e-014},
              {1.141056828884990e-009, -9.289835742028532e-010},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {2.228285989630589e+000, 1.814142639751731e+000}
      };
      const double Coeff_m7[size_Coeff_m][2] = {
              {-8.929405559006038e-011, -8.931137480031157e-011},
              {6.319814700961324e-005, 6.320244393309693e-005},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {4.785131212048377e-004, -5.877524860249395e-004},
              {4.271922830906078e+000, -5.247218808668549e+000}
      };
      const double Coeff_v7[size_Coeff_v][2] = {
              {-1.418731402291282e-012, 1.576782807097003e-011},
              {-5.512224505288543e-006, 1.914006058179041e-004},
              {1.000000000000000e+000, 1.000000000000000e+000},
              {-5.638714069888806e-006, 1.959753272178233e-004},
              {1.006201804172733e+000, 1.087101027065273e+000}
      };

//         *mixture = Calloc(1, GMRFLib_lgamma_mixture_tp);
//         (*mixture)->weights[0] = 0;
//         (*mixture)->means[0] = 0;
//         (*mixture)->variances[0] = 0;
//         (*mixture)->ncomp = 0;

         /*
          * Check input
          */
//         GMRFLib_ASSERT(n > 0, GMRFLib_EPARAMETER);
        if(n<1 && n>30000)
        {
        	cout<< " n is out side the range [1,30000] " << endl;
        }
        int nr = (int) n;
//
//         if (n < range[2][0] && nr != n) {
//                 GMRFLib_ASSERT(n < range[2][0] && (int) n != n, GMRFLib_EPARAMETER);
//                 return GMRFLib_EPARAMETER;
//         }

         /*
          * Mean and standard deviation of log-gamma
          */
         double mu = - gsl_sf_psi((double) n);
         double sigma = sqrt(gsl_sf_psi_1((double) n));

//         cout << "mu " << mu <<endl;
//         cout << "sigma " << sigma <<endl;
         /*
          * Single component
          */
         if (n > range[6][1]) {
                 vWeights[0] = 1.0;
                 vMeans[0] = mu;
                 vVariances[0] = sigma*sigma;
                 iNumOfComp[0] = 1;

                 return 1;
         }

         /*
          * Find appropriate range
          */
         int rtake = 0, ir, ic, nc, jc;

         for (ir = 0; ir < nrange; ir++) {
                 if (range[ir][0] <= n && n < range[ir][1] + 1) {
                         rtake = ir;
                 }
         }

         /*
          * Set pointers to coefficients
          */
         const double *P = NULL, *M = NULL, *V = NULL, *ptr = NULL;

         nr = nr - range[rtake][0];
         switch (rtake) {
         case 0:
                 P = &P10[nr][0];
                 M = &M10[nr][0];
                 V = &V10[nr][0];
                 break;
         case 1:
                 P = &P9[nr][0];
                 M = &M9[nr][0];
                 V = &V9[nr][0];
                 break;
         case 2:
                 P = &Coeff_p3[0][0];
                 M = &Coeff_m3[0][0];
                 V = &Coeff_v3[0][0];
                 break;
         case 3:
                 P = &Coeff_p4[0][0];
                 M = &Coeff_m4[0][0];
                 V = &Coeff_v4[0][0];
                 break;
         case 4:
                 P = &Coeff_p5[0][0];
                 M = &Coeff_m5[0][0];
                 V = &Coeff_v5[0][0];
                 break;
         case 5:
                 P = &Coeff_p6[0][0];
                 M = &Coeff_m6[0][0];
                 V = &Coeff_v6[0][0];
                 break;
         case 6:
                 P = &Coeff_p7[0][0];
                 M = &Coeff_m7[0][0];
                 V = &Coeff_v7[0][0];
                 break;
         default:
        	 	 cout << "Couldnt find interval" << endl;
                 return 0;
         }

         /*
          * Store number of components
          */
         nc = ncomp[rtake];
         iNumOfComp[0] = nc;

         /*
          * Explicit mixture parameters
          */
         if (rtake <= 1) {
                 for (ic = 0; ic < nc; ic++) {
                         vWeights[ic] = *P++;
                         vMeans[ic] = (*M++) * sigma + mu;
                         vVariances[ic] = (*V++) * sigma*sigma;
                 }

                 return 1;
         }

         /*
          * Mixture parameters by rational approximation
          */
         double numer = 0, sum;

         for (ic = 0; ic < nc; ic++) {
                 /*
                  * Weights
                  */
                 ptr = P++;
                 for (sum = jc = 0; jc < size_Coeff_p; jc++) {
                         sum = sum * n + (*ptr);
                         if (*ptr == 1) {
                                 numer = sum;
                                 sum = 0;
                         }
                         ptr = ptr + nc;
                 }
                 vWeights[ic] = numer / sum;

                 /*
                  * Means
                  */
                 ptr = M++;
                 for (sum = jc = 0; jc < size_Coeff_m; jc++) {
                         sum = sum * n + (*ptr);
                         if (*ptr == 1) {
                                 numer = sum;
                                 sum = 0;
                         }

                         ptr = ptr + nc;
                 }
                 vMeans[ic] = numer / sum * sigma + mu;

                 /*
                  * Variances
                  */
                 ptr = V++;
                 for (sum = jc = 0; jc < size_Coeff_v; jc++) {
                         sum = sum * n + (*ptr);
                         if (*ptr == 1) {
                                 numer = sum;
                                 sum = 0;
                         }

                         ptr = ptr + nc;
                 }
                 vVariances[ic]  = numer / sum * sigma*sigma;


         }

 #undef size_Coeff_p
 #undef size_Coeff_m
 #undef size_Coeff_v
 #undef size_Coeff_m3
 #undef size_Coeff_v3

         return 1;
 }

double LogitTransform(double x )
{
	return log(x/(1-x));
}

double LogitTransformBack(double x )
{
	return exp(x)/(1+exp(x));
}

double GammaLogPosterior(const gsl_vector *v, void *ParamPointer)
{
	double dValue=0;
	struct ObjectiveFunParamGamma * Param = (struct ObjectiveFunParamGamma *)ParamPointer;
	double dLogGamma= gsl_vector_get(v, 0);
	double dGamma=LogitTransformBack(dLogGamma);

//	cout << "XXXXXXXXXXXXXXXX Gamma iNumOF obs "<<  (Param->iNumOfObs)[0] << endl;
//	for(int i=0; i<(Param->iNumOfObs)[0]; i++)
//	{
				//cout << "vData in post" << i <<" is  "<< vData[i]<< endl;
				//cout << "vIndicator in post " << i <<" is  "<< (Param->vIndicator)[i]<< endl;
				//cout << "vSkellam in post " << i <<" is  "<< (Param->vSkellam)[i]<< endl;
//	}

	 for(int i=0; i<(Param->iNumOfObs)[0]; i++)
	 {
		 dValue+=log( dGamma*((Param->vIndicator)[i]) +(1-dGamma)*((Param->vSkellam)[i]) );
	 }
	 //cout <<"Function: " << dValue << " at: " << dGamma << endl;
	 //cout << "Prior A " << (Param->dPriorA)[0] << " Prior B " << (Param->dPriorB)[0] << endl;
	 dValue+=((Param->dPriorA)[0] -1)*log(dGamma)+((Param->dPriorB)[0] -1)*log(1-dGamma)  ;
	 return  -dValue/((Param->iNumOfObs)[0]);
}

/* The gradient of f, df = (df/dx, df/dy). */
void GammaLogPosteriorDerivative(const gsl_vector *v, void *ParamPointer,
       gsl_vector *df)
{
	double dValue=0;
	struct ObjectiveFunParamGamma * Param = (struct ObjectiveFunParamGamma *)ParamPointer;
	double dLogGamma= gsl_vector_get(v, 0);
	double dGamma=LogitTransformBack(dLogGamma);

	 for(int i=0; i<(Param->iNumOfObs)[0]; i++)
	 {
		 dValue+=((Param->vIndicator)[i]-(Param->vSkellam)[i] )/( dGamma*((Param->vIndicator)[i]) +(1-dGamma)*((Param->vSkellam)[i]) );
	 }

	 dValue+=((Param->dPriorA)[0] -1)/dGamma-((Param->dPriorB)[0] -1)/(1-dGamma)  ;

	 //cout <<"Derivative: " << -dValue <<" at: "<< dGamma << endl;
	 gsl_vector_set(df, 0, -dValue/((Param->iNumOfObs)[0]));

}

/* Compute both f and df together. */
void GammaLogPosteriorANDDerivative (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
  *f = GammaLogPosterior(x, params);
  GammaLogPosteriorDerivative(x, params, df);
}

void DrawGammaLaplaceApprox(  int iNumOfObs, const gsl_rng * gsl_random_num, struct AllParam sParam, double* vIndicator, double* vSkellam )
{
	/* Setting up the optimizer */
	int iter = 0;
	int max_iter=100;
	int status;

//	cout << "Num Of obs : " << iNumOfObs << endl;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf F;

	struct ObjectiveFunParamGamma sObjectiveFunParam;

	/* Passing parameters to the objective function */
	sObjectiveFunParam.dPriorA=sParam.dPriorGammaA ;
	sObjectiveFunParam.dPriorB= sParam.dPriorGammaB;
	sObjectiveFunParam.iNumOfObs=&iNumOfObs;
	sObjectiveFunParam.vIndicator=vIndicator;
	sObjectiveFunParam.vSkellam=vSkellam;

//	cout << "Num Of obs 2: " << iNumOfObs << endl;
	/* setting up the objective function */
	F.n = 1;
	F.f =&GammaLogPosterior;
	F.df = &GammaLogPosteriorDerivative;
	F.fdf = &GammaLogPosteriorANDDerivative;
	F.params = &sObjectiveFunParam;

	/* Starting point*/
	x = gsl_vector_alloc (1);
	gsl_vector_set (x,0, LogitTransform(0.15));

	/* Setting up the minimizer */
	T =  gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, 1);

	gsl_multimin_fdfminimizer_set (s, &F, x, 0.000001, 0.000001);

	/* Iterating the optimizer */
	do
	{
	      iter++;
	      status = gsl_multimin_fdfminimizer_iterate (s);

//	      cout<<iter<<endl;
	      if (status)
	        break;

	      status = gsl_multimin_test_gradient (s->gradient, 1e-2);

//	      if (status == GSL_SUCCESS)
//	      printf ("Minimum found at:\n");
//
//	      printf ("%5d %.5f %10.5f\n", iter,
//	              gsl_vector_get (s->x, 0),
//	              s->f);

	  }
	  while (status == GSL_CONTINUE && iter < max_iter);


//	  /* Mean */
	  double dMuLogit = gsl_vector_get (s->x, 0);


	  gsl_vector *dfx;
	  gsl_vector *dfxh;
	  gsl_vector *xh;
	  dfx=gsl_vector_alloc (1);
	  dfxh=gsl_vector_alloc (1);
	  xh=gsl_vector_alloc (1);
	  double dTol=0.0000001;

	  //cout <<  "Logit dMu:  "<<gsl_vector_get (s->x, 0) << endl;
	  //cout << "dMu: "<<LogitTransformBack(gsl_vector_get (s->x, 0)) << endl;

	  //cout <<  "Logit dMu +h: "<<gsl_vector_get (s->x, 0)+dTol << endl;
	  //cout << "dMu +h: "<<LogitTransformBack(gsl_vector_get (s->x, 0)+dTol) << endl;

	  gsl_vector_set (xh,0, gsl_vector_get (s->x, 0)+dTol);


	  GammaLogPosteriorDerivative(s->x, &sObjectiveFunParam ,dfx);
	  GammaLogPosteriorDerivative(xh, &sObjectiveFunParam ,dfxh);

	  /* Sigma */

	  double dSigmaLogit=sqrt(dTol/((gsl_vector_get (dfxh, 0)-gsl_vector_get (dfx, 0))*(iNumOfObs)));

	  //cout<< "dfxh: " << gsl_vector_get (dfxh,0) << endl;
	  //cout<< "dfx: " << gsl_vector_get (dfx,0)<< endl;
	  //cout<< "dSigmaLog: " << dSigmaLogit << endl;
	  //cout<< "dMu: " << dMuLogit << endl;
	  /* Draw new Gamma log */
	  //double dNewGammaLogit=dMuLogit+gsl_ran_gaussian(gsl_random_num,dSigmaLogit);
	  double dNewGammaLogit=dMuLogit+gsl_ran_tdist(gsl_random_num,4)*dSigmaLogit;

	  gsl_vector * newx;
	  newx=gsl_vector_alloc(1);
	  gsl_vector_set (newx,0, dNewGammaLogit);

	  gsl_vector * oldx;
	  oldx=gsl_vector_alloc(1);
	  gsl_vector_set (oldx,0, LogitTransform(sParam.dGamma[0]));

	  //cout<< "new Gamma logit: " << dNewGammaLogit << endl;
//	  cout<< "new Gamma: " << LogitTransformBack(dNewGammaLogit) << "  old Gamma: " << LogitTransformBack(gsl_vector_get (oldx, 0))<< endl;


	  /*Calculate Acceptance rate */
	  double dLogU=log(gsl_rng_uniform(gsl_random_num));
	  //cout << "Num Of obs 3 : " << iNumOfObs << endl;
//	  cout<<"Logpost new: "<< -GammaLogPosterior(newx, &sObjectiveFunParam )*iNumOfObs  <<
//	       " Logpost old: "<< +GammaLogPosterior(oldx, &sObjectiveFunParam )*iNumOfObs <<
//		 " Proposal old: "<< log(gsl_ran_gaussian_pdf(LogitTransform(sParam.dGamma[0])-dMuLogit,dSigmaLogit)) <<
//			 " Proposal new: "<< -log(gsl_ran_gaussian_pdf(dNewGammaLogit-dMuLogit,dSigmaLogit)) << endl;

//	 cout<<"Accept rate: "<<-GammaLogPosterior(newx, &sObjectiveFunParam )*iNumOfObs+ log(gsl_ran_gaussian_pdf(LogitTransform(sParam.dGamma[0])-dMuLogit,dSigmaLogit))
//							 +GammaLogPosterior(oldx, &sObjectiveFunParam)*iNumOfObs- log(gsl_ran_gaussian_pdf(dNewGammaLogit-dMuLogit,dSigmaLogit)) <<endl;


	  //double dAccept=fmin(-GammaLogPosterior(newx, &sObjectiveFunParam )*iNumOfObs+ log(gsl_ran_gaussian_pdf(LogitTransform(sParam.dGamma[0])-dMuLogit,dSigmaLogit))
		//			 +GammaLogPosterior(oldx, &sObjectiveFunParam)*iNumOfObs- log(gsl_ran_gaussian_pdf(dNewGammaLogit-dMuLogit,dSigmaLogit)),0);

	  double dNewLog=-GammaLogPosterior(newx, &sObjectiveFunParam )*iNumOfObs+ log(gsl_ran_tdist_pdf((LogitTransform(sParam.dGamma[0])-dMuLogit)/dSigmaLogit,4));
	  double dOldLog=+GammaLogPosterior(oldx, &sObjectiveFunParam)*iNumOfObs- log(gsl_ran_tdist_pdf((dNewGammaLogit-dMuLogit)/dSigmaLogit,4));
	  double dAccept=fmin(dNewLog+dOldLog,0);

	  if((dLogU<=dAccept ) & (dNewLog==dNewLog))
	  {
		  sParam.dGamma[0]=LogitTransformBack(dNewGammaLogit);

	  }
	  //cout  <<" dLogU " << dLogU<<" dAccept " << dAccept <<endl;






	  //cout << sParam.dGamma[0] << endl;
	  gsl_multimin_fdfminimizer_free (s);
	  gsl_vector_free (x);





}

void DrawGammaAdaptiveRW(  int iNumOfObs, int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam, double* vIndicator, double* vSkellam, double* dCovar,  double* dSum  )
{


	/* Proposal */
	double dOmega1=0.05;
	double dU=gsl_rng_uniform(gsl_random_num);
	double dNewGammaLogit;
	double dSigma;
	if(iNumOfIter<=1)
	{
		dSigma=0.1;
	}
	else if(iNumOfIter==2)
	{
		dSigma=pow(dCovar[0]-0.5*dSum[0]*dSum[0],0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}
	else
	{
		dSigma=pow(dCovar[0], 0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}

	if(dU<=dOmega1)
	{
		/* Static proposal */
		dNewGammaLogit=LogitTransform(sParam.dGamma[0])+gsl_ran_gaussian(gsl_random_num,0.1);
	}
	else
	{
		/* Proposal based on empirical covariance */
		dNewGammaLogit=LogitTransform(sParam.dGamma[0])+gsl_ran_gaussian(gsl_random_num,dSigma*2.38);
	}


	/* Acceptance Rate */
	gsl_vector * newx;
	newx=gsl_vector_alloc(1);
	gsl_vector_set (newx,0, dNewGammaLogit);

	gsl_vector * oldx;
	oldx=gsl_vector_alloc(1);
	gsl_vector_set (oldx,0, LogitTransform(sParam.dGamma[0]));


	/*Calculate Acceptance rate */
	struct ObjectiveFunParamGamma sObjectiveFunParam;

		/* Passing parameters to the objective function */
	sObjectiveFunParam.dPriorA=sParam.dPriorGammaA ;
	sObjectiveFunParam.dPriorB= sParam.dPriorGammaB;
	sObjectiveFunParam.iNumOfObs=&iNumOfObs;
	sObjectiveFunParam.vIndicator=vIndicator;
	sObjectiveFunParam.vSkellam=vSkellam;
	double dLogU=log(gsl_rng_uniform(gsl_random_num));
    double dAccept=fmin(-GammaLogPosterior(newx, &sObjectiveFunParam )*iNumOfObs+GammaLogPosterior(oldx, &sObjectiveFunParam)*iNumOfObs,0);

	if((dLogU<=dAccept) & (isthisfinite(LogitTransformBack(dNewGammaLogit))))
	{
		sParam.dGamma[0]=LogitTransformBack(dNewGammaLogit);
	}



	/* Covariance */
	if(iNumOfIter<=1)
	{
		dCovar[0]=dCovar[0]+pow(LogitTransform(sParam.dGamma[0]),2);
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);

	}
	else if(iNumOfIter==2)
	{
		dCovar[0]=(iNumOfIter-1)*pow(dSigma,2)/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));
	}
	else
	{
/* 	dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1)); */
/*		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+(dSum[0]/iNumOfIter)*(dSum[0]/iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-(dSum[0]/iNumOfIter)*(dSum[0]/(iNumOfIter+1)); */
		
		double dSum_old = dSum[0];
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0] = dCovar[0] + (LogitTransform(sParam.dGamma[0]) - dSum_old/iNumOfIter)*(LogitTransform(sParam.dGamma[0]) - dSum[0]/(iNumOfIter+1));
 		dCovar[0] = dCovar[0]/iNumOfIter;		
	}

	cout << "dSum " << dSum[0]<< endl;
	cout << "dCovar " << dCovar[0]<< endl;
	cout << "proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
	cout << "accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;


}

void MatrixTrans(double * mA,   int iRows,  int iCols, double * mC)
{

	for(int i=0; i<iRows; i++)
	{
		for(int j=0; j<iCols; j++)
		{
			mC[j*iRows+i]=mA[i*iCols+j];

		}
	}

}

void MatrixAdd(double * mA, double * mB,  int iRows,  int iCols, double * mC)
{
	for(int i=0; i<iRows; i++)
	{
		for(int j=0; j<iCols; j++)
		{
			mC[i*iCols+j]=mA[i*iCols+j]+mB[i*iCols+j];

		}
	}

}

void MatrixScalarMulti(double * mA, double dScalar, int iRows,  int iCols, double * mC)
{
	for(int i=0; i<iRows; i++)
	{
		for(int j=0; j<iCols; j++)
		{
			mC[i*iCols+j]=dScalar*mA[i*iCols+j];

		}
	}

}

void MatrixSub(double * mA, double * mB, int iRows,  int iCols, double * mC)
{
	for(int i=0; i<iRows; i++)
	{
		for(int j=0; j<iCols; j++)
		{
			mC[i*iCols+j]=mA[i*iCols+j]-mB[i*iCols+j];

		}
	}

}

void MatrixMulti(double * mA, double * mB, int iARows,  int iACols,   int iBRows,  int iBCols, double * mC)
{
	for(int i=0; i<iARows; i++)
	{
		for(int j=0; j<iBCols; j++)
		{
			mC[i*iBCols+j]=0;
			for(int k=0; k<iACols; k++)
			{
				//cout << " i " << i <<" j " << j << endl;
				//cout<< " k " << k << " mA[i*iARows+k] " << mA[i*iACols+k] <<" mB[k*iBCols+j] " << mB[k*iBCols+j]<<endl;
				mC[i*iBCols+j]+=mA[i*iACols+k]*mB[k*iBCols+j];
				//cout << " mC[i*iBCols+j] " <<mC[i*iBCols+j] << endl;
			}

		}
	}

}

void MatirxInv2x2(double* mA , int iSize, double * mC)
{

	if(iSize==1)
	{
		mC[0]=1/mA[0];
	}
	else if(iSize==2)
	{
		double dDet=mA[0]*mA[3]-mA[1]*mA[2];
		mC[0]=mA[3]/dDet;
		mC[1]=-mA[1]/dDet;
		mC[2]=-mA[2]/dDet;
		mC[3]=mA[0]/dDet;
	}
	else
	{
		cout<<"Error in inverse calculation"<< endl;
	}



}

void MatrixSandwitch(double * mOut, double * mIn, int iOutRows,int iOutCols,   int iInRows, int iInCols, double * mC )
{
	double * mTemp = new double[iInRows*iOutCols];
	double * mOutTrans=new double[iOutRows*iOutCols];
	MatrixTrans( mOut, iOutRows,iOutCols, mOutTrans);

//	cout<<"iOutRows "<<iOutRows<< endl;
//	cout<<"iOutCols "<<iOutCols<< endl;
//	cout<<"iInRows "<<iInRows<< endl;
//	cout<<"iInCols "<<iInCols<< endl;
//
//	cout<<"mOutTrans[0] "<<mOutTrans[0]<< endl;
//	cout<<"mOutTrans[1] "<<mOutTrans[1]<< endl;
//	cout<<"mOutTrans[2] "<<mOutTrans[2]<< endl;
//	cout<<"mOutTrans[3] "<<mOutTrans[3]<< endl;
//
//
//	cout<<"mOut[0] "<<mOut[0]<< endl;
//	cout<<"mOut[1] "<<mOut[1]<< endl;
//	cout<<"mOut[2] "<<mOut[2]<< endl;
//	cout<<"mOut[3] "<<mOut[3]<< endl;

	MatrixMulti(mIn, mOut,iInRows,iInCols, iOutRows,  iOutCols, mTemp);

//	cout<<"mTemp[0] "<<mTemp[0]<< endl;
//	cout<<"mTemp[1] "<<mTemp[1]<< endl;
//	cout<<"mTemp[2] "<<mTemp[2]<< endl;
//	cout<<"mTemp[3] "<<mTemp[3]<< endl;


	MatrixMulti(mOutTrans, mTemp, iOutCols, iOutRows, iInRows, iOutCols, mC);

	delete [] mTemp;
	delete [] mOutTrans;

}

void SystemMatrices(struct AllParam sParam, double* mZ, double* mT, double *mQ, double* mA1, double* mP1)
{
	/* Z */
	mZ[0]=1;
	mZ[1]=1;
	mZ[2]=1;
	mZ[3]=1;

//	for(int i=0; i<4; i++)
//	{
//		cout<<"mZ at "<<i<<" is "<<mZ[i]<<endl;
//	}

	/* T */
	mT[0]=1;
	mT[1]=0;
	mT[2]=0;
	mT[3]=sParam.dPhi[0];

//	for(int i=0; i<4; i++)
//	{
//		cout<<"mT at "<<i<<" is "<<mT[i]<<endl;
//	}

	/* Q */
	mQ[0]=0;
	mQ[1]=0;
	mQ[2]=0;
	mQ[3]=sParam.dSigma2[0];

//	for(int i=0; i<4; i++)
//	{
//		cout<<"mQ at "<<i<<" is "<<mQ[i]<<endl;
//	}

	/* A1 */
	mA1[0]=sParam.dPriorMuMean[0];
	mA1[1]=0;

//	for(int i=0; i<2; i++)
//	{
//		cout<<"mA1 at "<<i<<" is "<<mA1[i]<<endl;
//	}

	/* P1 */
	mP1[0]=sParam.dPriorMuSigma2[0];
	mP1[1]=0;
	mP1[2]=0;
	mP1[3]=sParam.dSigma2[0]/(1-sParam.dPhi[0]*sParam.dPhi[0]);

//	for(int i=0; i<4; i++)
//	{
//		cout<<"mP1 at "<<i<<" is "<<mP1[i]<<endl;
//	}




}

void KalmanFilter(double* mY, int* mN, double* mH, double* mZ, double* mT, double *mQ, double* mA1, double* mP1,  int iDimOfObs,
	 int iDimOfStates,  int iNumOfObs, double* mA, double * mP, double* mFInv, double* mV, double* mL)
{

	double * mZTrans= new double[iDimOfObs*iDimOfStates];


	for(int d=0; d<iDimOfStates; d++)
	{
		mA[d]=mA1[d];
		for(int k=0; k<iDimOfStates; k++)
		{
			mP[d*iDimOfStates+k]=mP1[d*iDimOfStates+k];
//			cout<<"mP[0] "<< d*iDimOfStates+k <<" " <<mP[d*iDimOfStates+k]<<endl;
		}
	}

	double* mNextA;
	double* mNextP;



	double * mZa= new double[iDimOfObs];
	double * mZPZ= new double[iDimOfObs*iDimOfObs];
	double * mK= new double[iDimOfObs*iDimOfStates];
	double * mPZFInv = new double[iDimOfStates*iDimOfObs];
	double * mZFInv = new double[iDimOfStates*iDimOfObs];
	double * mF= new double[iDimOfObs*iDimOfObs];
	double * mKZ= new double[iDimOfObs*iDimOfStates];
	double * mTa= new double[iDimOfStates];
	double * mKv= new double[iDimOfStates];
	double * mLTrans=new double[iDimOfStates*iDimOfStates];
	double * mPL=new double[iDimOfStates*iDimOfStates];
	double * mTPL=new double[iDimOfStates*iDimOfStates];

	for(int i=0; i<iNumOfObs; i++)
	{


		/* Setting the observation dimension to 1 when there is only one observation */
		if(mN[i]==0)
		{
			iDimOfObs=1;
		}

		MatrixTrans( mZ, iDimOfObs, iDimOfStates, mZTrans);

		/* v_t */
		if(i==0)
		{
			MatrixMulti(mZ, mA1,iDimOfObs,iDimOfStates, iDimOfStates, 1, mZa);
		}
		else
		{
			MatrixMulti(mZ, mA,iDimOfObs,iDimOfStates, iDimOfStates, 1, mZa);
		}

		MatrixSub(mY, mZa, iDimOfObs, 1, mV);

//		cout<<"XXXXXXXXX"<< "NumOfObs " << i<<" XXXXXXXX"<<endl;
//		cout<<"XXXXXXXXX"<< "DimOfObs " << i<<" XXXXXXXX"<<endl;
//		cout<<"mV[0] "<<mV[0]<<endl;
//		cout<<"mV[1] "<<mV[1]<<endl;
//
//		cout<<"mZ[0] "<<mZ[0]<< endl;
//		cout<<"mZ[1] "<<mZ[1]<< endl;
//		cout<<"mZ[2] "<<mZ[2]<< endl;
//		cout<<"mZ[3] "<<mZ[3]<< endl;
//
//		cout<<"mP[0] "<<mP[0]<< endl;
//		cout<<"mP[1] "<<mP[1]<< endl;
//		cout<<"mP[2] "<<mP[2]<< endl;
//		cout<<"mP[3] "<<mP[3]<< endl;
//
//		cout<<"mZTrans[0] "<<mZTrans[0]<< endl;
//		cout<<"mZTrans[1] "<<mZTrans[1]<< endl;
//		cout<<"mZTrans[2] "<<mZTrans[2]<< endl;
//		cout<<"mZTrans[3] "<<mZTrans[3]<< endl;
//
//		cout<<"mA[0] "<<mA[0]<<endl;
//		cout<<"mA[1] "<<mA[1]<<endl;
//
//		cout<<"mY[0] "<<mY[0]<<endl;
//		cout<<"mY[1] "<<mY[1]<<endl;
//


		/* F_t */
		if(i==0)
		{
			MatrixSandwitch(mZTrans, mP1,iDimOfStates,iDimOfObs, iDimOfStates, iDimOfStates, mZPZ);
		}
		else
		{
			MatrixSandwitch(mZTrans, mP,iDimOfStates,iDimOfObs, iDimOfStates, iDimOfStates, mZPZ);
		}

		MatrixAdd(mZPZ, mH, iDimOfObs, iDimOfObs, mF);

//		cout<<"mZ[0] "<<mZ[0]<< endl;
//		cout<<"mZ[1] "<<mZ[1]<< endl;
//		cout<<"mZ[2] "<<mZ[2]<< endl;
//		cout<<"mZ[3] "<<mZ[3]<< endl;

//		cout<<"mZPZ[0] "<<mZPZ[0]<< endl;
//		cout<<"mZPZ[1] "<<mZPZ[1]<< endl;
//		cout<<"mZPZ[2] "<<mZPZ[2]<< endl;
//		cout<<"mZPZ[3] "<<mZPZ[3]<< endl;

		/* K_t */
		MatirxInv2x2(mF , iDimOfObs,mFInv);
		MatrixMulti(mZTrans, mFInv,iDimOfStates, iDimOfObs, iDimOfObs, iDimOfObs, mZFInv);
		MatrixMulti(mP, mZFInv,iDimOfStates, iDimOfStates, iDimOfStates, iDimOfObs, mPZFInv);
		MatrixMulti(mT, mPZFInv,iDimOfStates, iDimOfStates, iDimOfStates, iDimOfObs, mK);

//		cout<<"mF[0] "<<mF[0]<< endl;
//		cout<<"mF[1] "<<mF[1]<< endl;
//		cout<<"mF[2] "<<mF[2]<< endl;
//		cout<<"mF[3] "<<mF[3]<< endl;
//
//		cout<<"mFInv[0] "<<mFInv[0]<< endl;
//		cout<<"mFInv[1] "<<mFInv[1]<< endl;
//		cout<<"mFInv[2] "<<mFInv[2]<< endl;
//		cout<<"mFInv[3] "<<mFInv[3]<< endl;

		/* L_t */
		MatrixMulti(mK, mZ,iDimOfStates, iDimOfObs, iDimOfObs, iDimOfStates, mKZ);
		MatrixSub(mT, mKZ, iDimOfStates, iDimOfStates, mL);

//		cout<<"mK[0] "<<mK[0]<< endl;
//		cout<<"mK[1] "<<mK[1]<< endl;
//		cout<<"mK[2] "<<mK[2]<< endl;
//		cout<<"mK[3] "<<mK[3]<< endl;



		if(i<iNumOfObs-1)
		{
			mNextA=mA+iDimOfStates;
			mNextP=mP+iDimOfStates*iDimOfStates;
			/* A_t */
			MatrixMulti(mT, mA,iDimOfStates, iDimOfStates, iDimOfStates, 1, mTa);
			MatrixMulti(mK, mV,iDimOfStates, iDimOfObs, iDimOfObs, 1, mKv);
			MatrixAdd(mTa, mKv, iDimOfStates, 1, mNextA);

			/* P_t */
			MatrixTrans( mL, iDimOfStates, iDimOfStates, mLTrans);
			MatrixMulti(mP, mLTrans,iDimOfStates, iDimOfStates, iDimOfStates,  iDimOfStates, mPL);
			MatrixMulti(mT, mPL,iDimOfStates, iDimOfStates, iDimOfStates,  iDimOfStates, mTPL);
			MatrixAdd(mTPL, mQ, iDimOfStates, iDimOfStates, mNextP);

//			cout<<"mA[0] "<<mA[0]<<endl;
//			cout<<"mA[1] "<<mA[1]<<endl;


			/* Set Dim back */
			if(mN[i]==0)
			{
				iDimOfObs=2;
			}

			/* Update pointers */
			mA=mA+iDimOfStates;
			mP=mP+iDimOfStates*iDimOfStates;
			mV=mV+iDimOfObs;
			mFInv=mFInv+iDimOfObs*iDimOfObs;
			mL=mL+iDimOfStates*iDimOfStates;
			mH=mH+iDimOfObs*iDimOfObs;
			mY=mY+iDimOfObs;
		}

//		cout<<"XXXXXXXXXXXXXXXXX"<<endl;

	}

	delete [] mZTrans;
	delete [] mZa;
	delete [] mZPZ;
	delete [] mK;
	delete [] mPZFInv;
	delete [] mZFInv;
	delete [] mF;
	delete [] mKZ;
	delete [] mTa;
	delete [] mKv;
	delete [] mLTrans;
	delete [] mPL;
	delete [] mTPL;
}

void KalmanSmoother(int * vN, double* mZ,  double * mA, double *mP, double *mV, double* mFInv, double* mL, int iDimOfObs,
		int iDimOfStates, unsigned int iNumOfObs, double *mAHat, double* mVHat )
{
	double * mZTrans= new double[iDimOfObs*iDimOfStates];


	double * mR= new double[iDimOfStates];
	double * mNextR= new double[iDimOfStates];
	double * mN= new double[iDimOfStates*iDimOfStates];
	double * mNextN = new double[iDimOfStates*iDimOfStates];

	double * mLTrans=new double[iDimOfStates*iDimOfStates];
	double * mFv =new double[iDimOfObs];
	double * mZFv=new double[iDimOfStates];
	double * mLR=new double[iDimOfStates];
	double * mPR=new double[iDimOfStates];
	double * mZFZ= new double[iDimOfStates*iDimOfStates];
	double * mLNL= new double[iDimOfStates*iDimOfStates];
	double * mPNP= new double[iDimOfStates*iDimOfStates];

	for(int i=0; i<iDimOfStates; i++)
	{
		mR[i]=0;
		for(int j=0; j<iDimOfStates; j++)
		{
			mN[i*iDimOfStates+j]=0;
		}
	}

	/* Set the pointer to the end of the arrays (backward recursion)*/
	mP=mP+(iNumOfObs-1)*iDimOfStates*iDimOfStates;
	mAHat=mAHat+(iNumOfObs-1)*iDimOfStates;
	mVHat=mVHat+(iNumOfObs-1)*iDimOfStates*iDimOfStates;
	mA=mA+(iNumOfObs-1)*iDimOfStates;
	mV=mV+(iNumOfObs-1)*iDimOfObs;
	mFInv=mFInv+(iNumOfObs-1)*iDimOfStates*iDimOfStates;
	mL=mL+(iNumOfObs-1)*iDimOfStates*iDimOfStates;

	for(int i=iNumOfObs-1; i>-1; i--)
	{
		if(vN[i]==0)
		{
			iDimOfObs=1;
		}

		MatrixTrans( mZ, iDimOfObs, iDimOfStates, mZTrans);
		/* mR */
		MatrixTrans( mL, iDimOfStates, iDimOfStates, mLTrans);
		MatrixMulti(mFInv, mV,iDimOfObs, iDimOfObs, iDimOfObs, 1, mFv);
		MatrixMulti(mZTrans, mFv,iDimOfStates, iDimOfObs, iDimOfObs, 1, mZFv);
		MatrixMulti(mLTrans, mR,iDimOfStates, iDimOfStates, iDimOfStates,1, mLR);
		MatrixAdd( mZFv, mLR, iDimOfStates, 1, mNextR);

		/* AHat */
		MatrixMulti(mP, mNextR,iDimOfStates, iDimOfStates, iDimOfStates,1, mPR);
		MatrixAdd( mA, mPR, iDimOfStates, 1, mAHat);

		/* mN*/
		MatrixSandwitch(mZ, mFInv,iDimOfObs,iDimOfStates, iDimOfObs, iDimOfObs, mZFZ);
		MatrixSandwitch(mL, mN,iDimOfStates,iDimOfStates, iDimOfStates, iDimOfStates, mLNL);
		MatrixAdd( mZFZ, mLNL, iDimOfStates, iDimOfStates, mNextN);

		/* VHat */
		MatrixSandwitch(mP, mNextN,iDimOfStates,iDimOfStates, iDimOfStates, iDimOfStates, mPNP);
		MatrixSub( mP, mPNP, iDimOfStates, iDimOfStates, mVHat);


		if(vN[i]==0)
		{
			iDimOfObs=2;
		}


		/* Update R and N values */
		for(int k=0; k<iDimOfStates; k++)
		{
			mR[k]=mNextR[k];
			for(int j=0; j<iDimOfStates; j++)
			{
				mN[k*iDimOfStates+j]=mNextN[k*iDimOfStates+j];
			}
		}
		if(i>0)
		{
			/* Update pointers */
			mP=mP-iDimOfStates*iDimOfStates;
			mAHat=mAHat-iDimOfStates;
			mVHat=mVHat-iDimOfStates*iDimOfStates;
			mA=mA-iDimOfStates;
			mV=mV-iDimOfObs;
			mFInv=mFInv-iDimOfStates*iDimOfStates;
			mL=mL-iDimOfStates*iDimOfStates;
		}

	}

	delete [] mZTrans;
	delete [] mR;
	delete [] mNextR;
	delete [] mN;
	delete [] mNextN ;
	delete [] mLTrans;
	delete [] mFv ;
	delete [] mZFv;
	delete [] mLR;
	delete [] mPR;
	delete [] mZFZ;
	delete [] mLNL;
	delete [] mPNP;

}

void Determinant2x2(double* mA, unsigned int iSize, double* dDet)
{
	if(iSize==1)
	{
		dDet[0]=mA[0];
	}
	else if(iSize==2)
	{
		dDet[0]=mA[0]*mA[3]-mA[1]*mA[2];
	}
	else
	{
		cout<<"Error in determinant calculation"<< endl;
	}
}

void CalculateLL(struct AllParam sParam, int iNumOfObs,  int iDimOfObs,
		int iDimOfStates, double * dLL )
{





	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatrices(sParam, mZ, mT, mQ,  mA1,  mP1);

	double * mA =new double[iNumOfObs*iDimOfStates];
	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
	double * mV =new double[iNumOfObs*iDimOfObs];
	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];
	KalmanFilter(sParam.mAuxY, sParam.vN, sParam.mAuxH,   mZ, mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mY at " << i<< " is "<< sParam.mAuxY[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*4; i++)
//	{
//		cout <<"mH at " << i<< " is "<< sParam.mAuxH[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//			cout <<"mV at " << i<< " is "<< mV[i]<< endl;
//	}

//	string sYFile="mY.csv";
//	WriteOutDoubleArray(sParam.mAuxY, iNumOfObs, 2, sYFile);
//	string sHFile="mH.csv";
//	WriteOutDoubleArray(sParam.mAuxH, iNumOfObs, 4, sHFile);
//	string sNFile="mN.csv";
//	WriteOutIntArray(sParam.vN, iNumOfObs, 1, sNFile);



	dLL[0]=0;
//	cout << "LL at init is " << dLL[0] << endl;
	double * dDetFInv= new double;
	double * dvFv=new double;


	for(int i=0; i<iNumOfObs; i++)
	{
		if(sParam.vN[i]==0)
		{
			iDimOfObs=1;
		}

		Determinant2x2(mFInv,iDimOfObs, dDetFInv);
		MatrixSandwitch(mV, mFInv,iDimOfObs,1, iDimOfObs, iDimOfObs, dvFv);
		dLL[0]=dLL[0]-0.5*iDimOfObs*log(2*PI)+0.5*(log(dDetFInv[0])-dvFv[0]);
//		cout << "LL Cont at " << i << " is " << 0.5*(log(dDetFInv[0])-dvFv[0]) << endl;
//		cout << "LL at " << i << " is " << dLL[0] << endl;

		if(sParam.vN[i]==0)
		{
			iDimOfObs=2;
		}
		if(i<iNumOfObs-1)
		{
			mV=mV+iDimOfObs;
			mFInv=mFInv+iDimOfObs*iDimOfObs;
		}
	}


	delete [] mZ;
	delete [] mT;
	delete [] mQ;
	delete [] mA1;
	delete [] mP1 ;
	delete [] mA;
	delete [] mP ;
	mFInv=mFInv-(iNumOfObs-1)*iDimOfObs*iDimOfObs;
	mV=mV-(iNumOfObs-1)*iDimOfObs;
	delete [] mFInv;
	delete [] mV;
	delete [] mL;
	delete dDetFInv;
	delete dvFv;


}

int Cholesky2x2(double* mA,  int iSize, double* mL)
{
	double dDet=mA[0]*mA[3]-mA[1]*mA[2];

	if(iSize==1)
	{
		if(mA[0]!=0)
		{
			mL[0]=sqrt(mA[0] );

			return 1;
		}
		else
		{
			cout<<"Error in cholesky"<< endl;
			cout<<"mA[0] "<< mA[0] <<endl;
			return 0;
		}


	}
	else if(iSize==2)
	{
		if(dDet>0)
		{
			mL[0]=sqrt(mA[0] );
			mL[1]=0;
			mL[2]=mA[2]/mL[0];
			mL[3]=sqrt(mA[3]-pow(mL[2],2));

			return 1;
		}
		else
		{
			cout<<"Error in cholesky"<< endl;
			cout<<"mA[0] "<< mA[0] <<endl;
			cout<<"mA[1] "<< mA[1] <<endl;
			cout<<"mA[2] "<< mA[2] <<endl;
			cout<<"mA[3] "<< mA[3] <<endl;
			return 0;
		}
	}
	else
	{
		cout<<"Error in cholesky"<< endl;
		return 0;

	}
}

void SSFRecursion(struct AllParam sParam, int * vN, double *mH, double* mZ, double* mT, double *mQ, double* mA1, double* mP1,
		int iNumOfObs,int iDimOfObs,int iDimOfStates, const gsl_rng * gsl_random_num, double* mYPlus, double* mAPlus)
{



	double * mQChol =new double[iDimOfStates*iDimOfStates];
	double * mHChol =new double[iDimOfObs*iDimOfObs];
	double * mP1Chol =new double[iDimOfStates*iDimOfStates];

	mQChol[0]=0;
	mQChol[1]=0;
	mQChol[2]=0;
	mQChol[3]=sqrt(sParam.dSigma2[0]);

//	Cholesky2x2(mP1,  iDimOfStates, mP1Chol);

	mP1Chol[0]=sqrt(sParam.dPriorMuSigma2[0]);
	mP1Chol[1]=0;
	mP1Chol[2]=0;
	mP1Chol[3]=sqrt(sParam.dSigma2[0]/(1-sParam.dPhi[0]*sParam.dPhi[0]));


//	cout<<"mQ[0] "<<mQ[0]<<endl;
//	cout<<"mQ[1] "<<mQ[1]<<endl;
//	cout<<"mQ[2] "<<mQ[2]<<endl;
//	cout<<"mQ[3] "<<mQ[3]<<endl;




	double * mEps =new double[iDimOfObs];
	double * mEta =new double[iDimOfStates];
	double * mEpsDraw =new double[iDimOfObs];
	double * mEtaDraw =new double[iDimOfStates];
	double * mTa =new double[iDimOfStates];
	double * mNextAPlus =new double[iDimOfStates];
	double * mZa =new double[iDimOfObs];
	double * mNextYPlus =new double[iDimOfObs];

	/* Random values to test SSF */
//	double * mRandom=new double[4*(iNumOfObs+1)];
//	mRandom[2]=0;
//	mRandom[3]=0;
//	mRandom[4*iNumOfObs]=0;
//	mRandom[4*iNumOfObs+1]=0;


	for(int i=0; i<iNumOfObs; i++)
	{
		if(vN[i]==0)
		{
			iDimOfObs=1;
			//mHChol[3]=0;
		}
//		cout << "SSF Recuriosn H Chol START " <<endl;
//		cout<<"vN[i] "<<vN[i]<<endl;

		Cholesky2x2(mH, iDimOfObs, mHChol);

//		cout<<"mH[0] "<<mH[0]<<endl;
//		cout<<"mH[3] "<<mH[3]<<endl;
//		cout<<"mHChol[0] "<<mHChol[0]<<endl;
//		cout<<"mHChol[1] "<<mHChol[1]<<endl;
//		cout<<"mHChol[2] "<<mHChol[2]<<endl;
//		cout<<"mHChol[3] "<<mHChol[3]<<endl;

//		cout << "SSF Recuriosn H Chol END " <<endl;
		/*Draw random numbers */
		for (int k=0;k<2;k++)
		{
			mEpsDraw[k]=gsl_ran_gaussian(gsl_random_num,1);

		}
		for (int k=0;k<2;k++)
		{
			mEtaDraw[k]=gsl_ran_gaussian(gsl_random_num,1);

		}




		/* Recursion */
		if(i==0)
		{
			MatrixMulti(mP1Chol, mEtaDraw,iDimOfStates, iDimOfStates, iDimOfStates, 1, mEta);
			//cout << "mEta[0] " <<mEta[0] <<endl;
			MatrixAdd(mA1, mEta, iDimOfStates, 1, mNextAPlus);
			//cout << "mAPlus[0] " <<mNextAPlus[0] <<endl;
		}
		else
		{
//			cout<<"mQChol[0] "<<mQChol[0]<<endl;
//			cout<<"mQChol[1] "<<mQChol[1]<<endl;
//			cout<<"mQChol[2] "<<mQChol[2]<<endl;
//			cout<<"mQChol[3] "<<mQChol[3]<<endl;


			MatrixMulti(mQChol, mEtaDraw,iDimOfStates, iDimOfStates, iDimOfStates, 1, mEta);

//			cout<<"Obs "<<i <<" mEtaDraw[0] "<<mEtaDraw[0]<<endl;
//			cout<<"Obs "<<i <<" mEtaDraw[1] "<<mEtaDraw[1]<<endl;
//			cout<<"Obs "<<i <<" mEta[0] "<<mEta[0]<<endl;
//			cout<<"Obs "<<i <<" mEta[1] "<<mEta[1]<<endl;

			MatrixMulti(mT, mAPlus,iDimOfStates, iDimOfStates, iDimOfStates, 1,mTa);
			MatrixAdd(mTa, mEta, iDimOfStates, 1, mNextAPlus);


		}

//		cout<<iDimOfObs<<endl;
//		cout<<"mH[0] "<<mH[0]<<endl;
//		cout<<"mH[1] "<<mH[1]<<endl;
//		cout<<"mH[2] "<<mH[2]<<endl;
//		cout<<"mH[3] "<<mH[3]<<endl;
//
//		cout<<"mHChol[0] "<<mHChol[0]<<endl;
//		cout<<"mHChol[1] "<<mHChol[1]<<endl;
//		cout<<"mHChol[2] "<<mHChol[2]<<endl;
//		cout<<"mHChol[3] "<<mHChol[3]<<endl;

		MatrixMulti(mHChol, mEpsDraw,iDimOfObs, iDimOfObs, iDimOfObs, 1, mEps);
		MatrixMulti(mZ, mNextAPlus,iDimOfObs, iDimOfStates, iDimOfStates, 1, mZa);
		MatrixAdd(mZa, mEps, iDimOfObs, 1, mNextYPlus);




		/* Delete Random later just for check */
//		mRandom[i*4]=mEta[0];
//		mRandom[i*4+1]=mEta[1];
//
//		if(iDimOfObs==1)
//		{
//			mRandom[(i+1)*4+2]=mEps[0];
//			mRandom[(i+1)*4+3]=0;
//		}
//		else
//		{
//			mRandom[(i+1)*4+2]=mEps[0];
//			mRandom[(i+1)*4+3]=mEps[1];
//		}

		if(vN[i]==0)
		{
			iDimOfObs=2;
		}

		/* Update pointers */
		if(i!=0 )
		{
			mAPlus=mAPlus+iDimOfStates;
			mYPlus=mYPlus+iDimOfObs;
		}

		if(i<iNumOfObs-1)
		{
			mH=mH+iDimOfObs*iDimOfObs;

		}
		for(int k=0; k<iDimOfStates; k++)
		{
			mAPlus[k]=mNextAPlus[k];

		}
		for(int k=0; k<iDimOfObs; k++)
		{
			mYPlus[k]=mNextYPlus[k];
		}
	}

	//cout << "mP1[0] " <<mP1[0] <<endl;

	delete [] mQChol ;
	delete [] mHChol  ;
	delete [] mP1Chol ;
	delete [] mEpsDraw;
	delete [] mEtaDraw;
	delete [] mEps;
	delete [] mEta;
	delete [] mZa;
	delete [] mTa;
	delete [] mNextAPlus;
	delete [] mNextYPlus;

	/* Save random values */
//	string sRFile="mRandom.csv";
//	WriteOutDoubleArray(mRandom, iNumOfObs+1, 4, sRFile);
//	cout<<"Done"<<endl;
//	delete [] mRandom;

}

void SimulationSmoother(struct AllParam sParam,double *mY,double *mH, int iNumOfObs,
		 int iDimOfObs,int iDimOfStates, const gsl_rng * gsl_random_num, double* mDraw)
{
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatrices(sParam,mZ, mT, mQ, mA1, mP1);

	double * mAHat =new double[iDimOfStates*iNumOfObs];
	double * mVHat =new double[iDimOfStates*iDimOfStates*iNumOfObs];
	double * mAHatPlus =new double[iDimOfStates*iNumOfObs];
	double * mVHatPlus =new double[iDimOfStates*iDimOfStates*iNumOfObs];
	double * mAPlus =new double[iDimOfStates*iNumOfObs];
	double * mYPlus =new double[iDimOfObs*iNumOfObs];

	double * mA =new double[iNumOfObs*iDimOfStates];
	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
	double * mV =new double[iNumOfObs*iDimOfObs];
	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];


	/* APlus */

//	string sYFile="mY.csv";
//	WriteOutDoubleArray(mY, iNumOfObs, 2, sYFile);
//	string sHFile="mH.csv";
//	WriteOutDoubleArray(mH, iNumOfObs, 4, sHFile);
//	string sNFile="mN.csv";
//	WriteOutIntArray(sParam.vN, iNumOfObs, 1, sNFile);


	SSFRecursion(sParam,sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1,  iNumOfObs,  iDimOfObs,iDimOfStates,  gsl_random_num, mYPlus, mAPlus);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAPlus at " << i<< " is "<< mAPlus[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mYPlus at " << i<< " is "<< mYPlus[i]<< endl;
//	}

	/* AHat */
	KalmanFilter(mY, sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
	KalmanSmoother( sParam.vN,mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHat,  mVHat );

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mY at " << i<< " is "<< mY[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mH at " << i<< " is "<< mH[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAHat at " << i<< " is "<< mAHat[i]<< endl;
//	}

	/* AHatPlus */

	KalmanFilter(mYPlus, sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mV at " << i<< " is "<< mV[i]<< endl;
//	}

	KalmanSmoother(sParam.vN, mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHatPlus,  mVHatPlus);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAHatPlus at " << i<< " is "<< mAHatPlus[i]<< endl;
//	}

	double * mAHatAHatPlus =new double[iDimOfStates*iNumOfObs];
	MatrixSub(mAHat, mAHatPlus, iNumOfObs, iDimOfStates, mAHatAHatPlus);
	MatrixAdd(mAHatAHatPlus, mAPlus, iNumOfObs, iDimOfStates, mDraw);

//	cout << "mDraw[0] " << mDraw[0] << endl;
//	cout << "mDraw[iNumOfObs-1] " << mDraw[2*iNumOfObs-1] << endl;

	delete [] mZ;
	delete [] mT;
	delete [] mQ;
	delete [] mA1;
	delete [] mP1;

	delete [] mAHat;
	delete [] mVHat;
	delete [] mAHatPlus;
	delete [] mVHatPlus;
	delete [] mAPlus;
	delete [] mYPlus;

	delete [] mA;
	delete [] mP;
	delete [] mFInv;
	delete [] mV;
	delete [] mL;
	delete [] mAHatAHatPlus;

}

void DrawXandMu(struct AllParam sParam, int iNumOfObs,const gsl_rng * gsl_random_num )
{
	double * mDraw=new double[2*iNumOfObs];

//	for(int i=0; i<iNumOfObs; i++)
//	{
//		sParam.mAuxY[2*i]=sParam.mAuxY[2*i]-sParam.vS[i];
//		sParam.mAuxY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS[i];
//	}



	SimulationSmoother(sParam,sParam.mAuxY,sParam.mAuxH, iNumOfObs,2,2, gsl_random_num, mDraw);

//	double * mY=new double[2*iNumOfObs];
//	double * mH=new double[4*iNumOfObs];
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		mY[2*i]=sParam.mAuxY[2*i]-sParam.vS[i];
//		mY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS[i];
//	}
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		mH[4*i]=sParam.mAuxH[4*i];
//		mH[4*i+1]=sParam.mAuxH[4*i+1];
//		mH[4*i+2]=sParam.mAuxH[4*i+2];
//		mH[4*i+3]=sParam.mAuxH[4*i+3];
//	}
//
//	SimulationSmoother(sParam,mY,sParam.mAuxH, iNumOfObs,2,2, gsl_random_num, mDraw);
//
//	delete [] mY;
//	delete [] mH;

	if((sParam.dMu[0]==sParam.dMu[0]) & (mDraw[0]!=mDraw[0]))
	{
		cout<<"Baj van :("<<endl;
		string sErrorYFile="vYError.csv";
		WriteOutDoubleArray( sParam.mAuxY , iNumOfObs, 1, sErrorYFile);

		sErrorYFile="vSError.csv";
		WriteOutDoubleArray( sParam.vS , iNumOfObs, 1, sErrorYFile);

		double * vParamError=new double[8];
		vParamError[0]=mDraw[0];
		vParamError[1]=sParam.dPhi[0];
		vParamError[2]=sParam.dSigma2[0];
		vParamError[3]=sParam.dGamma[0];
		vParamError[4]=sParam.vBeta[0];
		vParamError[5]=sParam.vBeta[1];
		vParamError[6]=sParam.vBeta[2];
		vParamError[7]=sParam.vBeta[3];

		string sErrorParamFile="vParamError.csv";
		WriteOutDoubleArray( vParamError , 1, 8, sErrorParamFile);


		delete [] vParamError;
	}

	sParam.dMu[0]=mDraw[0];

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vX[i]=mDraw[2*i+1];

	}

//	cout<<"mDraw[2*0] "<<mDraw[2*0]<<endl;
//	cout<<"mDraw[2*0+1] "<<mDraw[2*0+1]<<endl;
//	cout<<"mDraw[2*0+1] "<<mDraw[2*iNumOfObs-1]<<endl;
	delete [] mDraw;
}

double PhiSigmaLogPosteriorGSL(const gsl_vector *v, void *ParamPointer)
{

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double dLogitPhi= gsl_vector_get(v, 0);
	double dPhi=LogitTransformBack(dLogitPhi);
	double dLogSigma2= gsl_vector_get(v, 1);
	double dSigma2=exp(dLogSigma2);
	double dLLValue;

	double dPhiPrior=log(gsl_ran_beta_pdf ((dPhi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	double dSigmaPrior=log(gsl_ran_gamma_pdf (1/dSigma2, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));

//	cout <<"(Param->iNumOfObs)[0] in function GSL " <<  (Param->iNumOfObs)[0] << endl;
	/*LogLikelihood */
	double *dLL= new double;
	double dPhiTrue=(Param->dPhi)[0];
	double dSigma2True=(Param->dSigma2)[0];
	(Param->dPhi)[0]=dPhi;
	(Param->dSigma2)[0]=dSigma2;

	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLL );

	(Param->dSigma2)[0]=dSigma2True;
	(Param->dPhi)[0]=dPhiTrue;
	dLLValue=dLL[0];
	delete dLL;
//	cout <<"Value dPhi " <<dPhi << endl;
//	cout <<"Value dSigma " << dSigma2 << endl;
//	cout <<"Value " <<-(dLLValue+dPhiPrior+dSigmaPrior)/((Param->iNumOfObs)[0]) <<endl;
	return  -(dLLValue+dPhiPrior+dSigmaPrior)/((Param->iNumOfObs)[0]);


}

void PhiSigmaLogPosteriorDerivativeGSL(const gsl_vector *v, void *ParamPointer,
       gsl_vector *df)
{

	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double dLogitPhi= gsl_vector_get(v, 0);
	double dPhi=LogitTransformBack(dLogitPhi);
	double dHPhi=sqrt(dEps)*fmax( fabs(dLogitPhi),1);
	double xhPhi=dLogitPhi+dHPhi;
	double diffhPhi=xhPhi-dLogitPhi;
	double dPhiH=LogitTransformBack(xhPhi);
	double dLogSigma2= gsl_vector_get(v,1);
	double dSigma2=exp(dLogSigma2);
	double dHSigma=sqrt(dEps)*fmax(fabs(dLogSigma2),1);
	double xhSigma=dLogSigma2+dHSigma;
	double diffhSigma=xhSigma-dLogSigma2;
	double dSigmaH=exp(xhSigma);


//	cout <<"(Param->iNumOfObs)[0] in derivative GSL " <<  (Param->iNumOfObs)[0] << endl;
	/* Current LL*/
	/*LogLikelihood */
	double dPhiPrior=log(gsl_ran_beta_pdf ((dPhi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	double dSigmaPrior=log(gsl_ran_gamma_pdf (1/dSigma2, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	double *dLLCurrent= new double;

	double dPhiTrue=(Param->dPhi)[0];
	double dSigmaTrue=(Param->dSigma2)[0];
	(Param->dPhi)[0]=dPhi;
	(Param->dSigma2)[0]=dSigma2;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLLCurrent );

	cout<<"LL Current phi "<< (Param->dPhi)[0]<<" logit  "<< LogitTransform((Param->dPhi)[0])<<endl;
	cout<<"LL Current sigma "<< (Param->dSigma2)[0]<<" log "<< log((Param->dSigma2)[0])<<endl;
	double dAtCurrent=-(dLLCurrent[0]+dPhiPrior+dSigmaPrior)/((Param->iNumOfObs)[0]);

	delete dLLCurrent;


	/* Phi H LL*/
	/*LogLikelihood */
	double dPhiPriorH=log(gsl_ran_beta_pdf ((dPhiH+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	double *dLLPhiH= new double;

	(Param->dPhi)[0]=dPhiH;
	(Param->dSigma2)[0]=dSigma2;
	cout<<"LL Phi phi "<< (Param->dPhi)[0]<<" logit  "<< LogitTransform((Param->dPhi)[0])<<endl;
	cout<<"LL Phi sigma "<< (Param->dSigma2)[0]<<" log "<< log((Param->dSigma2)[0])<<endl;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLLPhiH );


	double dAtPhiH=-(dLLPhiH[0]+dPhiPriorH+dSigmaPrior)/((Param->iNumOfObs)[0]);

	delete dLLPhiH;

	/* Sigma H LL*/
	/*LogLikelihood */
	double dSigmaPriorH=log(gsl_ran_gamma_pdf (1/dSigmaH, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	double *dLLSigmaH= new double;

	(Param->dPhi)[0]=dPhi;
	(Param->dSigma2)[0]=dSigmaH;
	cout<<"LL Sigma phi "<< (Param->dPhi)[0]<<" logit  "<< LogitTransform((Param->dPhi)[0])<<endl;
	cout<<"LL Sigma sigma "<< (Param->dSigma2)[0]<<" log "<< log((Param->dSigma2)[0])<<endl;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLLSigmaH );


	double dAtSigmaH=-(dLLSigmaH[0]+dPhiPrior+dSigmaPriorH)/((Param->iNumOfObs)[0]);

	(Param->dPhi)[0]=dPhiTrue;
	(Param->dSigma2)[0]=dSigmaTrue;

	delete dLLSigmaH;

	cout <<"Derivative dPhi  " <<dPhi << endl;
	cout <<"Derivative dHPhi  " <<dHPhi << endl;
	cout <<"Derivative dPhi H  " <<dPhiH << endl;
	cout <<"Derivative dSigma " << dSigma2 << endl;
	cout <<"Derivative dHSigma " << dHSigma << endl;
	cout <<"Derivative dSigma H " << dSigmaH << endl;
	cout<<"Derivative AT LLCurrent "<< dAtCurrent <<endl ;
	cout<<"Derivative AT  LLPhi H "<< dAtPhiH <<endl ;
	cout<<"Derivative AT LLSigma H "<< dAtSigmaH <<endl ;
	cout<<"Derivative AT  DIFF Phi H "<<  dAtPhiH-dAtCurrent <<endl ;
	cout<<"Derivative AT DIFF Sigma H "<< dAtSigmaH-dAtCurrent <<endl ;
	cout<<"Derivative AT  Phi H "<<  (dAtPhiH-dAtCurrent)/diffhPhi <<endl ;
	cout<<"Derivative AT Sigma H "<< (dAtSigmaH-dAtCurrent)/diffhSigma<<endl ;
	gsl_vector_set(df, 0, (dAtPhiH-dAtCurrent)/diffhPhi);
	gsl_vector_set(df, 1,(dAtSigmaH-dAtCurrent)/diffhSigma);

}

void PhiSigmaLogPosteriorDerivativeGSL2(const gsl_vector *v, void *ParamPointer, gsl_vector *df)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= gsl_vector_get(v, 0);
	double y= gsl_vector_get(v, 1);

	double hforx=pow(dEps,1.0/3)*fmax( fabs(x),1);
	double hfory=pow(dEps,1.0/3)*fmax( fabs(y),1);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);


//	cout<<"DERIV phi0 "<< phi0 << endl;
//	cout<<"DERIV phi "<< phi << endl;
//	cout<<"DERIV phi1 "<< phi1 << endl;
//	cout<<"DERIV sigma0 "<< sigma0 << endl;
//	cout<<"DERIV sigma "<< sigma << endl;
//	cout<<"DERIV sigma1 "<< sigma1 << endl;


	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

	/* derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);


//	cout<<"DERIV fxh0 "<< fxh0 << endl;
//	cout<<"DERIV fxh1 "<< fxh1 << endl;
//	cout<<"DERIV hhforx "<< hhforx << endl;
//	cout<<"DERIV fxh1-fxh0 "<< fxh1-fxh0 << endl;
//	cout<<"DERIV (fxh1-fxh0)/hhforx "<< (fxh1-fxh0)/(2*hhforx) << endl;


	gsl_vector_set(df, 0, (fxh1-fxh0)/(2*hhforx));


	/* derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);


//	cout<<"DERIV fyh0 "<< fyh0 << endl;
//	cout<<"DERIV fyh1 "<< fyh1 << endl;
//	cout<<"DERIV hhfory "<< hhfory << endl;
//	cout<<"DERIV fyh1-fyh0 "<< fyh1-fyh0 << endl;
//	cout<<"DERIV (fyh1-fyh0)/hhfory "<< (fyh1-fyh0)/(2*hhfory) << endl;

	gsl_vector_set(df, 1,(fyh1-fyh0)/(2*hhfory));

	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;
//	cout<<"XXXXXXXXXXXXXXXXXXXXXX DERIV END XXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
}

void PhiSigmaLogPosteriorANDDerivativeGSL (const gsl_vector *v, void *ParamPointer,
        double *f, gsl_vector *df)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double dLogitPhi= gsl_vector_get(v, 0);
	double dPhi=LogitTransformBack(dLogitPhi);
	double dHPhi=sqrt(dEps)*fmax( fabs(dLogitPhi),1);
	double xhPhi=dLogitPhi+dHPhi;
	double diffhPhi=xhPhi-dLogitPhi;
	double dPhiH=LogitTransformBack(xhPhi);
	double dLogSigma= gsl_vector_get(v,1);
	double dSigma=exp(dLogSigma);
	double dHSigma=sqrt(dEps)*fmax(fabs(dLogSigma),1);
	double xhSigma=dLogSigma+dHSigma;
	double diffhSigma=xhSigma-dLogSigma;
	double dSigmaH=exp(xhSigma);


//	cout <<"(Param->iNumOfObs)[0] in derivative GSL " <<  (Param->iNumOfObs)[0] << endl;
	/* Current LL*/
	/*LogLikelihood */
	double dPhiPrior=log(gsl_ran_beta_pdf ((dPhi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	double dSigmaPrior=log(gsl_ran_gamma_pdf (1/dSigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	double *dLLCurrent= new double;

	double dPhiTrue=(Param->dPhi)[0];
	double dSigmaTrue=(Param->dSigma2)[0];
	(Param->dPhi)[0]=dPhi;
	(Param->dSigma2)[0]=dSigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLLCurrent );

	cout<<"LL Current phi "<< (Param->dPhi)[0]<<" logit  "<< LogitTransform((Param->dPhi)[0])<<endl;
	cout<<"LL Current sigma "<< (Param->dSigma2)[0]<<" log "<< log((Param->dSigma2)[0])<<endl;
	double dAtCurrent=-(dLLCurrent[0]+dPhiPrior+dSigmaPrior)/((Param->iNumOfObs)[0]);

	delete dLLCurrent;


	/* Phi H LL*/
	/*LogLikelihood */
	double dPhiPriorH=log(gsl_ran_beta_pdf ((dPhiH+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	double *dLLPhiH= new double;

	(Param->dPhi)[0]=dPhiH;
	(Param->dSigma2)[0]=dSigma;
	cout<<"LL Phi phi "<< (Param->dPhi)[0]<<" logit  "<< LogitTransform((Param->dPhi)[0])<<endl;
	cout<<"LL Phi sigma "<< (Param->dSigma2)[0]<<" log "<< log((Param->dSigma2)[0])<<endl;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLLPhiH );


	double dAtPhiH=-(dLLPhiH[0]+dPhiPriorH+dSigmaPrior)/((Param->iNumOfObs)[0]);

	delete dLLPhiH;

	/* Sigma H LL*/
	/*LogLikelihood */
	double dSigmaPriorH=log(gsl_ran_gamma_pdf (1/dSigmaH, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	double *dLLSigmaH= new double;

	(Param->dPhi)[0]=dPhi;
	(Param->dSigma2)[0]=dSigmaH;
	cout<<"LL Sigma phi "<< (Param->dPhi)[0]<<" logit  "<< LogitTransform((Param->dPhi)[0])<<endl;
	cout<<"LL Sigma sigma "<< (Param->dSigma2)[0]<<" log "<< log((Param->dSigma2)[0])<<endl;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  dLLSigmaH );


	double dAtSigmaH=-(dLLSigmaH[0]+dPhiPriorH+dSigmaPriorH)/((Param->iNumOfObs)[0]);

	(Param->dPhi)[0]=dPhiTrue;
	(Param->dSigma2)[0]=dSigmaTrue;

	delete dLLSigmaH;

	f[0]=dAtCurrent;
	gsl_vector_set(df, 0, (dAtPhiH-dAtCurrent)/diffhPhi);
	gsl_vector_set(df, 1,(dAtSigmaH-dAtCurrent)/diffhSigma);


}

void PhiSigmaLogPosteriorANDDerivativeGSL2 (const gsl_vector *v, void *ParamPointer,
        double *f, gsl_vector *df)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= gsl_vector_get(v, 0);
	double y= gsl_vector_get(v, 1);

	double hforx=pow(dEps,1.0/3)*fmax( fabs(x),1);
	double hfory=pow(dEps,1.0/3)*fmax( fabs(y),1);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);

	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

	/* function evaluation */
	double fx=0;
	fx=fx+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fx=fx+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fx=-(fx+LL[0])/((Param->iNumOfObs)[0]);

	f[0]=fx;

	/* derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);

	gsl_vector_set(df, 0, (fxh1-fxh0)/(2*hhforx));


	/* derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);

	gsl_vector_set(df, 1,(fyh1-fyh0)/(2*hhfory));

	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;

}

void CalculatingHessian(const gsl_vector *v, void *ParamPointer,double * mHessian)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= gsl_vector_get(v, 0);
	double y= gsl_vector_get(v, 1);


	double hforx=pow(dEps,1.0/4)*fmax( (double) fabs((double) x),1.0);
	double hfory=pow(dEps,1.0/4)*fmax( (double) fabs( (double) y),1.0);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);

	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

//	cout<<"HESSIAN phi0 "<< phi0 << endl;
//	cout<<"HESSIAN phi "<< phi << endl;
//	cout<<"HESSIAN phi1 "<< phi1 << endl;
//	cout<<"HESSIAN sigma0 "<< sigma0 << endl;
//	cout<<"HESSIAN sigma "<< sigma << endl;
//	cout<<"HESSIAN sigma1 "<< sigma1 << endl;

	/* function evaluation */
	double fx=0;
	fx=fx+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fx=fx+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fx=-(fx+LL[0])/((Param->iNumOfObs)[0]);

	/* 2nd derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);

	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);

//	cout<<"HESSIAN fxh0 "<< fxh0 << endl;
//	cout<<"HESSIAN fxh1 "<< fxh1 << endl;
//	cout<<"HESSIAN fx "<< fx << endl;
//	cout<<"HESSIAN hhforx "<< hhforx << endl;
//	cout<<"HESSIAN fxh1-2*fx+fxh0 "<< fxh1-2*fx+fxh0<< endl;
//	cout<<"HESSIAN (fxh1-2*fx+fxh0)/pow(hhforx,2)"<< (fxh1-2*fx+fxh0)/pow(hhforx,2)<< endl;

	mHessian[0]= (Param->iNumOfObs)[0]*((fxh1-2*fx+fxh0)/pow(hhforx,2));

	/* 2nd derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);

	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);

	mHessian[3]= (Param->iNumOfObs)[0]*((fyh1-2*fx+fyh0)/pow(hhfory,2));


	/* corss derivative with respect to x and y */
	/* x forward  y forward*/
	double fxh1yh1=0;
	fxh1yh1=fxh1yh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1yh1=fxh1yh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1yh1=-(fxh1yh1+LL[0])/((Param->iNumOfObs)[0]);
	/* x forward  y backward*/
	double fxh1yh0=0;
	fxh1yh0=fxh1yh0+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1yh0=fxh1yh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1yh0=-(fxh1yh0+LL[0])/((Param->iNumOfObs)[0]);
	/* x backward  y forward*/
	double fxh0yh1=0;
	fxh0yh1=fxh0yh1+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0yh1=fxh0yh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0yh1=-(fxh0yh1+LL[0])/((Param->iNumOfObs)[0]);
	/* x backward  y backward*/
	double fxh0yh0=0;
	fxh0yh0=fxh0yh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0yh0=fxh0yh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0yh0=-(fxh0yh0+LL[0])/((Param->iNumOfObs)[0]);

	mHessian[1]=(Param->iNumOfObs)[0]*((fxh1yh1-fxh1yh0-fxh0yh1+fxh0yh0)/(4*hhforx*hhfory));
	mHessian[2]=(Param->iNumOfObs)[0]*((fxh1yh1-fxh1yh0-fxh0yh1+fxh0yh0)/(4*hhforx*hhfory));


	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;


}

void CalculateLLVector(struct AllParam sParam, int iNumOfObs,  int iDimOfObs,
		int iDimOfStates, double * mLL )
{





	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatrices(sParam, mZ, mT, mQ,  mA1,  mP1);

	double * mA =new double[iNumOfObs*iDimOfStates];
	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
	double * mV =new double[iNumOfObs*iDimOfObs];
	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];
	KalmanFilter(sParam.mAuxY, sParam.vN, sParam.mAuxH,   mZ, mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mY at " << i<< " is "<< sParam.mAuxY[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*4; i++)
//	{
//		cout <<"mH at " << i<< " is "<< sParam.mAuxH[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//			cout <<"mV at " << i<< " is "<< mV[i]<< endl;
//	}

//	string sYFile="mY.csv";
//	WriteOutDoubleArray(sParam.mAuxY, iNumOfObs, 2, sYFile);
//	string sHFile="mH.csv";
//	WriteOutDoubleArray(sParam.mAuxH, iNumOfObs, 4, sHFile);
//	string sNFile="mN.csv";
//	WriteOutIntArray(sParam.vN, iNumOfObs, 1, sNFile);




//	cout << "LL at init is " << dLL[0] << endl;
	double * dDetFInv= new double;
	double * dvFv=new double;


	for(int i=0; i<iNumOfObs; i++)
	{
		if(sParam.vN[i]==0)
		{
			iDimOfObs=1;
		}

		Determinant2x2(mFInv,iDimOfObs, dDetFInv);
		MatrixSandwitch(mV, mFInv,iDimOfObs,1, iDimOfObs, iDimOfObs, dvFv);
		mLL[i]=-0.5*iDimOfObs*log(2*PI)+0.5*(log(dDetFInv[0])-dvFv[0]);
//		cout << "LL Cont at " << i << " is " << 0.5*(log(dDetFInv[0])-dvFv[0]) << endl;
//		cout << "LL at " << i << " is " << dLL[0] << endl;

		if(sParam.vN[i]==0)
		{
			iDimOfObs=2;
		}
		if(i<iNumOfObs-1)
		{
			mV=mV+iDimOfObs;
			mFInv=mFInv+iDimOfObs*iDimOfObs;
		}
	}


	delete [] mZ;
	delete [] mT;
	delete [] mQ;
	delete [] mA1;
	delete [] mP1 ;
	delete [] mA;
	delete [] mP ;
	mFInv=mFInv-(iNumOfObs-1)*iDimOfObs*iDimOfObs;
	mV=mV-(iNumOfObs-1)*iDimOfObs;
	delete [] mFInv;
	delete [] mV;
	delete [] mL;
	delete dDetFInv;
	delete dvFv;


}

void CalculatingHessianFormJacobian(const gsl_vector *v, void *ParamPointer,double * mHessian)
{
	double dEps=2.22044604925031308e-16;


		struct AllParam * Param = (struct AllParam *)ParamPointer;
		double x= gsl_vector_get(v, 0);
		double y= gsl_vector_get(v, 1);

		double hforx=pow(dEps,1.0/3)*fmax( fabs(x),1);
		double hfory=pow(dEps,1.0/3)*fmax( fabs(y),1);

		double xh1=x+hforx;
		double xh0=x-hforx;

		double yh1=y+hfory;
		double yh0=y-hfory;

		double hhforx=xh1-xh0;
		double hhfory=yh1-yh0;


		double phi1=LogitTransformBack(xh1);
		double phi0=LogitTransformBack(xh0);
		double phi=LogitTransformBack(x);
		double sigma1=exp(yh1);
		double sigma0=exp(yh0);
		double sigma=exp(y);


		cout<<"DERIV hforx "<< hforx<< endl;
		cout<<"DERIV hhforx "<< hhforx<< endl;
		cout<<"DERIV phi0 "<< phi0 << endl;
		cout<<"DERIV phi "<< phi << endl;
		cout<<"DERIV phi1 "<< phi1 << endl;
		cout<<"DERIV sigma0 "<< sigma0 << endl;
		cout<<"DERIV sigma "<< sigma << endl;
		cout<<"DERIV sigma1 "<< sigma1 << endl;


		double phiTrue=(Param->dPhi)[0];
		double sigmaTrue=(Param->dSigma2)[0];

		double* mLLforward= new double[(Param->iNumOfObs)[0]];
		double* mLLbackward= new double[(Param->iNumOfObs)[0]];
		double* mJacobian=new double[2*(Param->iNumOfObs)[0]];
		/* derivative with respect to x */
		/* forward */
		(Param->dPhi)[0]=phi1;
		(Param->dSigma2)[0]=sigma;
		CalculateLLVector( Param[0],(Param->iNumOfObs)[0], 2, 2,  mLLforward );
		/* backward */
		(Param->dPhi)[0]=phi0;
		(Param->dSigma2)[0]=sigma;
		CalculateLLVector( Param[0],(Param->iNumOfObs)[0], 2, 2, mLLbackward  );

		double dLLf=0;
		double dLLb=0;
		for(int i=0 ; i<(Param->iNumOfObs)[0]; i++)
		{
//			cout << "mLLforward at "<< i << " is " << mLLforward[i] << endl;
//			cout << "mLLbackward at "<< i << " is " << mLLbackward[i] << endl;
//			cout << "hhforx at "<< i << " is " << 2*hhforx << endl;
			dLLf=dLLf+mLLforward[i];
			dLLb=dLLb+mLLbackward[i];
			mJacobian[i*2]= ( mLLforward[i]- mLLbackward[i])/(2*hhforx);
//			cout << "mJacobian[i*2] at "<< i << " is " << mJacobian[i*2]<< endl;
		}


		cout <<"dLLf "<<dLLf <<endl;
		cout <<"dLLb "<<dLLb <<endl;
		/* derivative with respect to y */
		/* forward */
		(Param->dPhi)[0]=phi;
		(Param->dSigma2)[0]=sigma1;
		CalculateLLVector( Param[0],(Param->iNumOfObs)[0], 2, 2,  mLLforward );

		/* backward */
//		double fyh0=0;
		(Param->dPhi)[0]=phi;
		(Param->dSigma2)[0]=sigma0;
		CalculateLLVector( Param[0],(Param->iNumOfObs)[0], 2, 2, mLLbackward  );


		for(int i=0 ; i<(Param->iNumOfObs)[0]; i++)
		{
			mJacobian[i*2+1]= ( mLLforward[i]- mLLbackward[i])/(2*hhfory);
		}


		(Param->dPhi)[0]=phiTrue;
		(Param->dSigma2)[0]=sigmaTrue;
		delete [] mLLforward;
		delete [] mLLbackward;


		double * mJacobianTrans=new double[2*(Param->iNumOfObs)[0]];
		MatrixTrans( mJacobian, (Param->iNumOfObs)[0],2, mJacobianTrans);
		MatrixMulti(mJacobianTrans ,mJacobian,2, (Param->iNumOfObs)[0], (Param->iNumOfObs)[0], 2, mHessian);

		delete [] mJacobian;
		delete [] mJacobianTrans;

}

double BivariateStudentTDensity(double* mX, double* mMean, double* mScale, unsigned int iDf)
{
	double dDet= mScale[0]*mScale[3]-mScale[1]*mScale[2];
	double dPdf=tgamma((iDf+2.0)/2.0)/(sqrt(dDet*pow(PI*iDf,2))*tgamma(iDf/2.0));
	double * mDeMean=new double[2];
	MatrixSub(mX, mMean, 2, 1, mDeMean);
	double *mScaleInv=new double [4];
	MatirxInv2x2(mScale , 2,mScaleInv);
	double * dSumS=new double;
	MatrixSandwitch(mDeMean, mScaleInv,2,1, 2, 2, dSumS );
	dPdf*=pow(1.0+dSumS[0]/(double) iDf, -(iDf+2.0)/2.0);

	delete [] mDeMean;
	delete [] mScaleInv;
	delete dSumS;
	return dPdf;
}

void DrawPhiSigmaLaplaceApprox( struct AllParam sParam , int iNumOfObs, const gsl_rng * gsl_random_num )
{

	for(int i=0; i<iNumOfObs; i++)
	{

		if(!isthisfinite(sParam.mAuxY[2*i+1])| !isthisfinite( sParam.vS[i]) )
		{
			string sFile="vAuxYPhiError.csv";
			WriteOutDoubleArray( sParam.mAuxY , iNumOfObs, 2, sFile);

			sFile="vSPhiError.csv";
			WriteOutDoubleArray( sParam.vS , iNumOfObs, 1, sFile);
			exit(1);
		}

		sParam.mAuxY[2*i]=sParam.mAuxY[2*i]-sParam.vS[i];
		sParam.mAuxY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS[i];
	}


	/* Setting up the optimizer */
	int iter = 0;
	int max_iter=100;
	int status;

//	cout << "Num Of obs : " << iNumOfObs << endl;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf F;


	/* setting up the objective function */
	F.n = 2;
	F.f =&PhiSigmaLogPosteriorGSL;
	F.df = &PhiSigmaLogPosteriorDerivativeGSL2;
	F.fdf = &PhiSigmaLogPosteriorANDDerivativeGSL2;
	F.params = &sParam;

	/* Starting point*/
	x = gsl_vector_alloc (2);
	gsl_vector_set (x,0, LogitTransform(sParam.dPhi[0]));
	gsl_vector_set (x,1,log(sParam.dSigma2[0]));

	/* Setting up the minimizer */
	T =  gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, 2);

	gsl_multimin_fdfminimizer_set (s, &F, x, 0.000001, 0.000001);

	/* Iterating the optimizer */
	do
	{
	      iter++;
	      status = gsl_multimin_fdfminimizer_iterate (s);

//	      cout<<iter<<endl;
	      if (status)
	        break;

	      status = gsl_multimin_test_gradient (s->gradient, 1e-2);

//	      if (status == GSL_SUCCESS)
//	      printf ("Minimum found at:\n");
//
//	      printf ("%5d %.5f %10.5f\n", iter,
//	              gsl_vector_get (s->x, 0),
//	              s->f);

	  }
	  while (status == GSL_CONTINUE && iter < max_iter);


//	  /* Mean */
	  double* mMeanLogit=new double[2];
	  mMeanLogit[0]=gsl_vector_get (s->x, 0);
	  mMeanLogit[1]=gsl_vector_get(s->x,1);


	  double* mInfoMatrix=new double[4];
//	  double* mInfoMatrixJacobian=new double[4];
	  gsl_vector *atx;
	  atx=gsl_vector_alloc (2);
	  gsl_vector_set (atx,0, gsl_vector_get (s->x, 0));
	  gsl_vector_set (atx,1, gsl_vector_get (s->x, 1));
	  CalculatingHessian(atx, &sParam, mInfoMatrix);
//	  CalculatingHessianFormJacobian(atx, &sParam, mInfoMatrixJacobian);
//	  gsl_vector *deriv;
//	  deriv=gsl_vector_alloc (2);
//	  PhiSigmaLogPosteriorDerivativeGSL2(atx, &sParam, deriv);
//	  cout<<"derive 0 " <<  gsl_vector_get (deriv, 0) <<endl;
//	  cout<<"derive 1 " <<  gsl_vector_get (deriv, 1) <<endl;
//	  cout<<"Phi mean "<< LogitTransformBack(mMeanLogit[0])<<endl;
//	  cout<<"Sigma mean "<< exp(mMeanLogit[1])<<endl;




//	  gsl_vector *dfx;
//	  gsl_vector *dfxh;
//	  gsl_vector *xh;
//	  dfx=gsl_vector_alloc (2);
//	  dfxh=gsl_vector_alloc (2);
//	  xh=gsl_vector_alloc (2);
//	  double* mInfoMatrix=new double[4];
//
//	  cout<< "XXXXXXXXXXXX DERIVATIVES XXXXXXXXXXXXXXXXXXXXX"<<endl;
//	  cout<< "XXXXXXXXXXXX DERIVATIVES XXXXXXXXXXXXXXXXXXXXX"<<endl;
//	  cout<< "XXXXXXXXXXXX DERIVATIVES XXXXXXXXXXXXXXXXXXXXX"<<endl;
//
//	  double dEps=2.22044604925031308e-16;
//	  double dHPhi=sqrt(dEps)*fmax(fabs(mMeanLogit[0]),1);
//	  double xhPhi=mMeanLogit[0]+dHPhi;
//	  double diffhPhi=xhPhi-mMeanLogit[0];
//	  double dHSigma=sqrt(dEps)*fmax(fabs(mMeanLogit[1]),1);
//	  double xhSigma=mMeanLogit[1]+dHPhi;
//	  double diffhSigma=xhSigma-mMeanLogit[1];
//
//
//	  cout << "mMeanLogit[0] "<< mMeanLogit[0] <<endl;
//	  cout << "mMeanLogit[1] "<< mMeanLogit[1] <<endl;
//	  cout << "diffhPhi"<< diffhPhi <<endl;
//	  cout << "diffhSigma "<< diffhSigma <<endl;
//	  cout << "xhPhi"<< xhPhi <<endl;
//	  cout << "xhSigma "<< xhSigma <<endl;
//
//	  /* dPhi second derivative*/
//	  gsl_vector_set (xh,0, xhPhi);
//	  gsl_vector_set (xh,1, gsl_vector_get (s->x, 1));
//	  PhiSigmaLogPosteriorDerivativeGSL(s->x, &sParam ,dfx);
//	  cout << "d Phi Sigma 0 " << gsl_vector_get (dfx, 0) << endl;
//	  cout << "d Phi Sigma  1 " << gsl_vector_get (dfx, 1) << endl;
//
//	  PhiSigmaLogPosteriorDerivativeGSL(xh, &sParam ,dfxh);
//
//	  cout << "d Phi+h 0 " << gsl_vector_get (dfxh, 0) << endl;
//	  cout << "d Phi+h 1 " << gsl_vector_get (dfxh, 1) << endl;
//
//	  mInfoMatrix[0]=iNumOfObs*((gsl_vector_get (dfxh, 0)-gsl_vector_get (dfx, 0))/diffhPhi);
//	  mInfoMatrix[2]=iNumOfObs*((gsl_vector_get (dfxh, 1)-gsl_vector_get (dfx, 1))/diffhPhi);
//	  mInfoMatrix[1]=iNumOfObs*((gsl_vector_get (dfxh, 1)-gsl_vector_get (dfx, 1))/diffhPhi);
//	  /*dSigma second derivative */
//	  gsl_vector_set (xh,0, gsl_vector_get (s->x, 0));
//	  gsl_vector_set (xh,1, xhSigma);
//	  PhiSigmaLogPosteriorDerivativeGSL(xh, &sParam ,dfxh);
//	  cout << "d Sigma+h 0 " << gsl_vector_get (dfxh, 0) << endl;
//	  cout << "d Sigma+h 1 " << gsl_vector_get (dfxh, 1) << endl;
//
//	  mInfoMatrix[3]=iNumOfObs*((gsl_vector_get (dfxh, 1)-gsl_vector_get (dfx, 1))/diffhSigma);
	  /*dSigma dPhi derivative */
//	  mInfoMatrix[1]=iNumOfObs*((gsl_vector_get (dfxh, 0)-gsl_vector_get (dfx,0))/diffhSigma);
//	  mInfoMatrix[2]=iNumOfObs*((gsl_vector_get (dfxh, 0)-gsl_vector_get (dfx,0))/diffhSigma);

//	  for(int k=0; k<4; k++)
//	  {
//	  	   cout<<"mInfoMatrix at "<< k <<" is "<<mInfoMatrix[k]<<endl;
//	  	   cout<<"mInfoMatrixJacobian at "<< k <<" is "<<mInfoMatrixJacobian[k]<<endl;
//	  }

	  double dDet=mInfoMatrix[0]*mInfoMatrix[3]-mInfoMatrix[1]*mInfoMatrix[2];
	  if( (dDet<0 ) | ( dDet!=dDet))
	  {
		  mInfoMatrix[1]=0;
		  mInfoMatrix[2]=0;
	  }

	  if((mInfoMatrix[0]<0 )|( mInfoMatrix[0]!=mInfoMatrix[0]) )
	  {
		  mInfoMatrix[0]=0.01;
	  }
	  if( (mInfoMatrix[3]<0) |(mInfoMatrix[3]!=mInfoMatrix[3]) )
	  {
		  mInfoMatrix[3]=0.01;
	  }
	  /* Inverse info matrix */
	  double* mInvInfoMatrix=new double[4];
	  MatirxInv2x2(mInfoMatrix , 2,mInvInfoMatrix);
//	  if(dDet<0)
//	  {
//	 		  mInvInfoMatrix[1]=0.075;
//	 		  mInvInfoMatrix[2]=0.075;
//	  }
//	  for(int k=0; k<4; k++)
//	  {
//	  		  cout<<"mInvInfoMatrix at "<< k <<" is "<<mInvInfoMatrix[k]<<endl;
//	  }

	  /* Cholesky Inverse info matrix  */
	  double* mCholInvInfoMatrix=new double[4];
//	  cout<<"Phi Sigma Laplace chol start"<<endl;
	  Cholesky2x2(mInvInfoMatrix, 2, mCholInvInfoMatrix);
//	  cout<<"Phi Sigma Laplace chol end"<<endl;
//	  for(int k=0; k<4; k++)
//	  {
//	 	  		  cout<<"mCholInvInfoMatrix at "<< k <<" is "<<mCholInvInfoMatrix[k]<<endl;
//	  }


	  /* Draw new Phi and Sigma from t density */
	  double* mRandom=new double[2];
	  double* mNewLogit=new double[2];
	  unsigned int  iDf=5;
	  double dSqrtW=sqrt((double) iDf/gsl_ran_chisq (gsl_random_num, iDf));

	  for (int k=0;k<2;k++)
	  {
		  mRandom[k]=dSqrtW*gsl_ran_gaussian(gsl_random_num,1);

	  }
	  MatrixMulti(mCholInvInfoMatrix, mRandom,2, 2, 2, 1, mNewLogit);
	  MatrixAdd(mMeanLogit, mNewLogit, 2, 1, mNewLogit);



	  gsl_vector * newx;
	  newx=gsl_vector_alloc(2);
	  gsl_vector_set (newx,0, mNewLogit[0]);
	  gsl_vector_set (newx,1, mNewLogit[1]);

	  gsl_vector * oldx;
	  oldx=gsl_vector_alloc(2);
	  gsl_vector_set (oldx,0, LogitTransform(sParam.dPhi[0]));
	  gsl_vector_set (oldx,1, log(sParam.dSigma2[0]));
	  double * mOldLogit=new double[2];
	  mOldLogit[0]=LogitTransform(sParam.dPhi[0]);
	  mOldLogit[1]=log(sParam.dSigma2[0]);


	  /*Calculate Acceptance rate */
	  double dLogU=log(gsl_rng_uniform(gsl_random_num));

//	  double sigma_dPhi=sqrt(mInvInfoMatrix[0]);
//	  double sigma_dSigma=sqrt(mInvInfoMatrix[3]);
//	  double rho_dPhiSigma=mInvInfoMatrix[1]/(sigma_dPhi*sigma_dSigma);


//	  cout << "new Phi " <<LogitTransformBack( mNewLogit[0]) <<endl;
//	  cout << "new Sigma " <<exp( mNewLogit[1]) <<endl;
//	  cout << "LL at new "<<-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs <<endl;
//	  cout << "LL at old "<< +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs <<endl;
//	  cout << "proposal at new "<< log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)- mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma))<<endl;
//	  cout << "proposal at old "<< log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma)) <<endl;
//	  cout << "acceptance "<< -PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma))
//							 		 +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)- mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma)) <<endl;


//	  double dAccept=fmin(-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma))
//					 		 +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)-mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma)),0);

	  double dNewLog=-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(BivariateStudentTDensity(mOldLogit,  mMeanLogit, mInvInfoMatrix , iDf));
	  double dOldLog= -PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs+ log(BivariateStudentTDensity(mNewLogit,mMeanLogit, mInvInfoMatrix , iDf));
	  double dAccept=fmin(dNewLog-dOldLog,0);

//	  cout<< "g++"<<endl;
//	  cout << "new "<<dNewLog <<endl;
//	  cout << "old "<< dOldLog <<endl;
//	  cout  <<"dLogU " << dLogU<<" dAccept " << dAccept <<endl;



//	  cout <<"new phi "<<LogitTransformBack(mNewLogit[0])<<endl;
//	  cout <<"old phi "<<LogitTransformBack(mOldLogit[0])<<endl;
//	  cout <<"new phi log"<<-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(BivariateStudentTDensity(mOldLogit,  mMeanLogit, mInvInfoMatrix , iDf))<<endl;
//	  cout <<"old phi log "<<+PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(BivariateStudentTDensity(mNewLogit,mMeanLogit, mInvInfoMatrix , iDf))<<endl;
//	  cout <<"Accept "<<dAccept<<endl;
	  if((dLogU<=dAccept) & (LogitTransformBack(mNewLogit[0])!=1) & (dNewLog==dNewLog)  )
	  {

//		  printf ("new phi: %1.12f  \n", LogitTransformBack(mNewLogit[0]) );
		  sParam.dPhi[0]=LogitTransformBack(mNewLogit[0]);
		  sParam.dSigma2[0]=exp(mNewLogit[1]);

	  }
//	  cout  <<"Accepted Phi " << sParam.dPhi[0] <<" Sigma "<< sParam.dSigma2[0] << endl;






	  //cout << sParam.dGamma[0] << endl;
	  gsl_multimin_fdfminimizer_free (s);
	  gsl_vector_free (x);

	  delete [] mInfoMatrix;
	  delete [] mInvInfoMatrix;
	  delete [] mCholInvInfoMatrix;
	  delete [] mRandom;
	  delete [] mNewLogit;
	  delete [] mMeanLogit;



}

void CalculatingGradientSGD(double * mLogit, void *ParamPointer, double * mGradient)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= mLogit[0];
	double y= mLogit[1];

	double hforx=pow(dEps,1.0/3)*fmax( (double) fabs( (double) x), 1.0);
	double hfory=pow(dEps,1.0/3)*fmax( (double) fabs( (double) y), 1.0);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);


//	cout<<"DERIV phi0 "<< phi0 << endl;
//	cout<<"DERIV phi "<< phi << endl;
//	cout<<"DERIV phi1 "<< phi1 << endl;
//	cout<<"DERIV sigma0 "<< sigma0 << endl;
//	cout<<"DERIV sigma "<< sigma << endl;
//	cout<<"DERIV sigma1 "<< sigma1 << endl;


	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

	/* derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);


//	cout<<"DERIV fxh0 "<< fxh0 << endl;
//	cout<<"DERIV fxh1 "<< fxh1 << endl;
//	cout<<"DERIV hhforx "<< hhforx << endl;
//	cout<<"DERIV fxh1-fxh0 "<< fxh1-fxh0 << endl;
//	cout<<"DERIV (fxh1-fxh0)/hhforx "<< (fxh1-fxh0)/(2*hhforx) << endl;


	mGradient[0]= (fxh1-fxh0)/(2*hhforx);


	/* derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);


//	cout<<"DERIV fyh0 "<< fyh0 << endl;
//	cout<<"DERIV fyh1 "<< fyh1 << endl;
//	cout<<"DERIV hhfory "<< hhfory << endl;
//	cout<<"DERIV fyh1-fyh0 "<< fyh1-fyh0 << endl;
//	cout<<"DERIV (fyh1-fyh0)/hhfory "<< (fyh1-fyh0)/(2*hhfory) << endl;

	mGradient[1]= (fyh1-fyh0)/(2*hhfory);

	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;
//	cout<<"XXXXXXXXXXXXXXXXXXXXXX DERIV END XXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
}

void CalculatingHessianSGD(double * mLogit, void *ParamPointer,double * mHessian)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x=mLogit[0];
	double y=mLogit[1];


	double hforx=pow(dEps,1.0/4)*fmax( (double) fabs((double) x),1.0);
	double hfory=pow(dEps,1.0/4)*fmax( (double) fabs((double) x),1.0);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);

	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

//	cout<<"HESSIAN phi0 "<< phi0 << endl;
//	cout<<"HESSIAN phi "<< phi << endl;
//	cout<<"HESSIAN phi1 "<< phi1 << endl;
//	cout<<"HESSIAN sigma0 "<< sigma0 << endl;
//	cout<<"HESSIAN sigma "<< sigma << endl;
//	cout<<"HESSIAN sigma1 "<< sigma1 << endl;

	/* function evaluation */
	double fx=0;
	fx=fx+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fx=fx+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fx=-(fx+LL[0])/((Param->iNumOfObs)[0]);

	/* 2nd derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);

//	cout<<"HESSIAN fxh0 "<< fxh0 << endl;
//	cout<<"HESSIAN fxh1 "<< fxh1 << endl;
//	cout<<"HESSIAN fx "<< fx << endl;
//	cout<<"HESSIAN hhforx "<< hhforx << endl;
//	cout<<"HESSIAN fxh1-2*fx+fxh0 "<< fxh1-2*fx+fxh0<< endl;
//	cout<<"HESSIAN (fxh1-2*fx+fxh0)/pow(hhforx,2) "<< (fxh1-2*fx+fxh0)/pow(hhforx,2)<< endl;

	mHessian[0]= ((fxh1-2*fx+fxh0)/pow(hhforx,2));

	/* 2nd derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);

	mHessian[3]=((fyh1-2*fx+fyh0)/pow(hhfory,2));


	/* corss derivative with respect to x and y */
	/* x forward  y forward*/
	double fxh1yh1=0;
	fxh1yh1=fxh1yh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1yh1=fxh1yh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1yh1=-(fxh1yh1+LL[0])/((Param->iNumOfObs)[0]);
	/* x forward  y backward*/
	double fxh1yh0=0;
	fxh1yh0=fxh1yh0+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1yh0=fxh1yh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh1yh0=-(fxh1yh0+LL[0])/((Param->iNumOfObs)[0]);
	/* x backward  y forward*/
	double fxh0yh1=0;
	fxh0yh1=fxh0yh1+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0yh1=fxh0yh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma1;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0yh1=-(fxh0yh1+LL[0])/((Param->iNumOfObs)[0]);
	/* x backward  y backward*/
	double fxh0yh0=0;
	fxh0yh0=fxh0yh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0yh0=fxh0yh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma0;
	CalculateLL( Param[0],(Param->iNumOfObs)[0], 2, 2,  LL );
	fxh0yh0=-(fxh0yh0+LL[0])/((Param->iNumOfObs)[0]);

	mHessian[1]=((fxh1yh1-fxh1yh0-fxh0yh1+fxh0yh0)/(4*hhforx*hhfory));
	mHessian[2]=((fxh1yh1-fxh1yh0-fxh0yh1+fxh0yh0)/(4*hhforx*hhfory));


	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;


}

void DrawPhiSigmaLaplaceApproxSGD( struct AllParam sParam ,  int iNumOfObs, const gsl_rng * gsl_random_num )
{
	int iNumOfIter=1;

	double* mOldLogit=new double[2];
	double* mPrevLogit=new double[2];
	double* mMeanLogit=new double[2];
	double* mNewLogit=new double[2];
	double* mGradient=new double[2];
	double* mHessian=new double[4];
	double* mInvHessian=new double[4];

	mOldLogit[0]=LogitTransform(sParam.dPhi[0]);
	mOldLogit[1]=log(sParam.dSigma2[0]);

	mPrevLogit[0]=mOldLogit[0];
	mPrevLogit[1]=mOldLogit[1];

	for(int i=0; i<iNumOfIter;i++)
	{
		/* Gradient*/
		CalculatingGradientSGD(mPrevLogit, &sParam ,mGradient);
		/* Hessian */
		CalculatingHessianSGD(mPrevLogit, &sParam ,mHessian);
		double dDetH=mHessian[0]*mHessian[3]-mHessian[1]*mHessian[2];
		if(dDetH<0)
		{
			mHessian[1]=0;
			mHessian[2]=0;
		}
		MatirxInv2x2(mHessian , 2,mInvHessian);
		for(int k=0; k<4; k++)
		{
			cout<<"mInvHessin at "<< k <<" is "<<mInvHessian[k]<<endl;
		}


//		mInvHessian[0]=150;
//		mInvHessian[1]=30;
//		mInvHessian[2]=30;
//		mInvHessian[3]=20;


		/* Update */
		MatrixMulti(mInvHessian,mGradient ,2, 2, 2, 1, mMeanLogit);
		MatrixSub( mOldLogit,mMeanLogit, 2, 1, mMeanLogit);

		mPrevLogit[0]=mMeanLogit[0];
		mPrevLogit[1]=mMeanLogit[1];

	}

	double* mInfoMatrix=new double[4];
//	double* mInfoMatrixJacobian=new double[4];
	gsl_vector *atx;
	atx=gsl_vector_alloc (2);
    gsl_vector_set (atx,0,mMeanLogit[0]);
	gsl_vector_set (atx,1, mMeanLogit[1]);
    CalculatingHessian(atx, &sParam, mInfoMatrix);
//    CalculatingHessianFormJacobian(atx, &sParam, mInfoMatrixJacobian);
	cout<<"Phi mean "<< LogitTransformBack(mMeanLogit[0])<<endl;
	cout<<"Sigma mean "<< exp(mMeanLogit[1])<<endl;

	  for(int k=0; k<4; k++)
	  {
		  cout<<"mInfoMatrix at "<< k <<" is "<<mInfoMatrix[k]<<endl;
//		  cout<<"mInfoMatrixJacobian at "<< k <<" is "<<mInfoMatrixJacobian[k]<<endl;
	  }
	  double dDet=mInfoMatrix[0]*mInfoMatrix[3]-mInfoMatrix[1]*mInfoMatrix[2];
	  if(dDet<0)
	  {
		  mInfoMatrix[1]=0;
		  mInfoMatrix[2]=0;
	  }
	  /* Inverse info matrix */
	  double* mInvInfoMatrix=new double[4];
	  MatirxInv2x2(mInfoMatrix , 2,mInvInfoMatrix);
//	  if(dDet<0)
//	  {
//	 		  mInvInfoMatrix[1]=0.075;
//	 		  mInvInfoMatrix[2]=0.075;
//	  }

	  if(mInvInfoMatrix[0]<0 || mInvInfoMatrix[3]<0 )
	  {
		  mInvInfoMatrix[1]=0.02;
		  mInvInfoMatrix[2]=0.005;
		  mInvInfoMatrix[1]=0.005;
		  mInvInfoMatrix[2]=0.003;
	  }

	  for(int k=0; k<4; k++)
	  {
	  		  cout<<"mInvInfoMatrix at "<< k <<" is "<<mInvInfoMatrix[k]<<endl;
	  }


	  /* Cholesky Inverse info matrix  */
	  double* mCholInvInfoMatrix=new double[4];
	  Cholesky2x2(mInvInfoMatrix, 2, mCholInvInfoMatrix);



	  /* Draw new Phi and Sigma from t density */
	  double* mRandom=new double[2];
	  unsigned int  iDf=5;
	  double dSqrtW=sqrt((double) iDf/gsl_ran_chisq (gsl_random_num, iDf));

	  for (int k=0;k<2;k++)
	  {
		  mRandom[k]=dSqrtW*gsl_ran_gaussian(gsl_random_num,1);

	  }
	  MatrixMulti(mCholInvInfoMatrix, mRandom,2, 2, 2, 1, mNewLogit);
	  MatrixAdd(mMeanLogit, mNewLogit, 2, 1, mNewLogit);



	  gsl_vector * newx;
	  newx=gsl_vector_alloc(2);
	  gsl_vector_set (newx,0, mNewLogit[0]);
	  gsl_vector_set (newx,1, mNewLogit[1]);

	  gsl_vector * oldx;
	  oldx=gsl_vector_alloc(2);
	  gsl_vector_set (oldx,0, LogitTransform(sParam.dPhi[0]));
	  gsl_vector_set (oldx,1, log(sParam.dSigma2[0]));



	  /*Calculate Acceptance rate */
	  double dLogU=log(gsl_rng_uniform(gsl_random_num));

//	  double sigma_dPhi=sqrt(mInvInfoMatrix[0]);
//	  double sigma_dSigma=sqrt(mInvInfoMatrix[3]);
//	  double rho_dPhiSigma=mInvInfoMatrix[1]/(sigma_dPhi*sigma_dSigma);


//	  cout << "new Phi " <<LogitTransformBack( mNewLogit[0]) <<endl;
//	  cout << "new Sigma " <<exp( mNewLogit[1]) <<endl;
//	  cout << "LL at new "<<-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs <<endl;
//	  cout << "LL at old "<< +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs <<endl;
//	  cout << "proposal at new "<< log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)- mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma))<<endl;
//	  cout << "proposal at old "<< log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma)) <<endl;
//	  cout << "acceptance "<< -PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma))
//							 		 +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)- mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma)) <<endl;


//	  double dAccept=fmin(-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma))
//					 		 +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)-mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma)),0);

	  double dAccept=fmin(-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(BivariateStudentTDensity(mOldLogit,  mMeanLogit, mInvInfoMatrix , iDf))
	  					 		 +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(BivariateStudentTDensity(mNewLogit,mMeanLogit, mInvInfoMatrix , iDf)),0);

	  cout  <<"dLogU " << dLogU<<" dAccept " << dAccept <<endl;
	  if(dLogU<=dAccept)
	  {
		  sParam.dPhi[0]=LogitTransformBack(mNewLogit[0]);
		  sParam.dSigma2[0]=exp(mNewLogit[1]);
	  }
	  cout  <<"Accepted Phi " << sParam.dPhi[0] <<" Sigma "<< sParam.dSigma2[0] << endl;






	  //cout << sParam.dGamma[0] << endl;

	  delete [] mInfoMatrix;
	  delete [] mInvInfoMatrix;
	  delete [] mCholInvInfoMatrix;
	  delete [] mRandom;
	  delete [] mNewLogit;
	  delete [] mMeanLogit;
	  delete [] mOldLogit;
	  delete [] mPrevLogit;
	  delete [] mGradient;
	  delete [] mHessian;
	  delete [] mInvHessian;


}

void PhiSigmaLogPosterior(struct AllParam sParam, int iNumOfObs, int iDimOfObs,
		 int iDimOfStates, double * dLogPost)
{
	/*Log Priors */
	/* Phi Prior */
	dLogPost[0]=log(gsl_ran_beta_pdf ((sParam.dPhi[0]+1)/2, sParam.dPriorPhiA[0], sParam.dPriorPhiB[0]));

//	cout<<"Phi prior "<<dLogPost[0]<<endl;
	/* Sigma2 Prior */
	dLogPost[0]=dLogPost[0]+log(gsl_ran_gamma_pdf (1/sParam.dSigma2[0], sParam.dPriorSigmaA[0], sParam.dPriorSigmaB[0]));
//	cout<<"Sigma prior "<<dLogPost[0]<<endl;
	/*LogLikelihood */
	double *dLL= new double;
	CalculateLL( sParam,iNumOfObs, iDimOfObs, iDimOfStates,  dLL );

	dLogPost[0]=dLogPost[0]+dLL[0];
//	cout<<"dLL "<<dLL[0]<<endl;
//	cout<<"full posterior "<<dLogPost[0]<<endl;
	delete dLL;

}

void DrawPhiSigmaAdaptiveRW( struct AllParam sParam,int iNumOfObs,  int iNumOfIter,
		const gsl_rng * gsl_random_num, double* mCovar,  double* mSum )
{
	/* Proposal */
	double dOmega1=0.05;
	double dU=gsl_rng_uniform(gsl_random_num);
//	double dNewPhiLogit;
	double* mChol=new double[4];
	double* mSigma=new double[4];
	double* mSumSum=new double[4];
	double* mParamParam=new double[4];
	double *mSumTrans= new double[2];
	double *mLogParam=new double[2];
	double *mLogParamTrans=new double[2];
	double *mNewLogParam=new double[2];
	double *mRandom=new double[2];
	int iFlag=1;

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.mAuxY[2*i]=sParam.mAuxY[2*i]-sParam.vS[i];
		sParam.mAuxY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS[i];
	}


	if(iNumOfIter<=1)
	{
//		cout<<"first two iter"<<endl;
		mChol[0]=0.1/sqrt(2);
		mChol[1]=0;
		mChol[2]=0;
		mChol[3]=0.1/sqrt(2);
	}
	else if(iNumOfIter==2)
	{
//		cout<<"third iter"<<endl;
		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
		MatrixScalarMulti(mSumSum,0.5*0.5,2, 2, mSumSum);
		MatrixScalarMulti(mCovar,0.5,2, 2, mCovar);
		MatrixSub(mCovar,mSumSum, 2,2, mSigma);
		iFlag=Cholesky2x2(mSigma,2, mChol);
//		cout<<"mSigma[0] "<<mSigma[0]<<endl;
//		cout<<"mSigma[1] "<<mSigma[1]<<endl;
//		cout<<"mSigma[2] "<<mSigma[2]<<endl;
//		cout<<"mSigma[3] "<<mSigma[3]<<endl;



	}
	else
	{
//		cout<<"later iter"<<endl;
		iFlag=Cholesky2x2(mCovar,2, mChol);

	}

	mLogParam[0]=LogitTransform(sParam.dPhi[0]);
	mLogParam[1]=log(sParam.dSigma2[0]);

	if(dU<=dOmega1)
	{
		/* Static proposal */
		mChol[0]=0.1/sqrt(2);
		mChol[1]=0;
		mChol[2]=0;
		mChol[3]=0.1/sqrt(2);

		mRandom[0]=gsl_ran_gaussian(gsl_random_num,1);
		mRandom[1]=gsl_ran_gaussian(gsl_random_num,1);

//		cout<<"mRandom[0] "<<mRandom[0]<<endl;
//		cout<<"mRandom[1] "<<mRandom[1]<<endl;

		MatrixMulti(mChol, mRandom,2, 2, 2, 1, mNewLogParam);
		MatrixAdd(mNewLogParam, mLogParam, 2, 1, mNewLogParam);


	}
	else
	{

		mRandom[0]=gsl_ran_gaussian(gsl_random_num,1);
		mRandom[1]=gsl_ran_gaussian(gsl_random_num,1);
//		cout<<"mRandom[0] "<<mRandom[0]<<endl;
//		cout<<"mRandom[1] "<<mRandom[1]<<endl;
//		cout<<"mLogParam[0] "<<mLogParam[0]<<endl;
//		cout<<"mLogParam[1] "<<mLogParam[1]<<endl;

		if(iFlag==0 )
		{
			mChol[0]=0.1/sqrt(2);
			mChol[1]=0;
			mChol[2]=0;
			mChol[3]=0.1/sqrt(2);
		}
		else
		{
			MatrixScalarMulti(mChol,2.38/sqrt(2),2, 2, mChol);
		}
//		cout<<"iFlag "<<iFlag<<endl;
//		cout<<"mChol[0] "<<mChol[0]<<endl;
//		cout<<"mChol[1] "<<mChol[1]<<endl;
//		cout<<"mChol[2] "<<mChol[2]<<endl;
//		cout<<"mChol[3] "<<mChol[3]<<endl;


//		mChol[0] =0.001;
//		mChol[1] =0;
//		mChol[2] =0;
//		mChol[3] =0.0005;

		MatrixMulti(mChol, mRandom,2, 2, 2, 1, mNewLogParam);
		MatrixAdd(mNewLogParam, mLogParam, 2, 1, mNewLogParam);

	}



	/* Acceptance Rate */

	double dLogU=log(gsl_rng_uniform(gsl_random_num));

	double* dNewLogPosterior= new double;
	double* dOldLogPosterior= new double;


	sParam.dPhi[0]=LogitTransformBack(mNewLogParam[0]);
	sParam.dSigma2[0]=exp(mNewLogParam[1]);

//	cout<<"XXXXX new XXXXXXXX"<<endl;
//	cout<<" new Phi " << sParam.dPhi[0]<< endl;
//	cout<<" new Sigma " << sParam.dSigma[0]<< endl;
	PhiSigmaLogPosterior(sParam, iNumOfObs, 2, 2,dNewLogPosterior);

	sParam.dPhi[0]=LogitTransformBack(mLogParam[0]);
	sParam.dSigma2[0]=exp(mLogParam[1]);

//	cout<<"XXXXX old XXXXXXXX"<<endl;
//	cout<<" old Phi " << sParam.dPhi[0]<< endl;
//	cout<<" old Sigma2 " << sParam.dSigma2[0]<< endl;
	PhiSigmaLogPosterior(sParam,iNumOfObs, 2, 2,dOldLogPosterior);

    double dAccept=fmin(dNewLogPosterior[0]-dOldLogPosterior[0],0);


//    cout<< "old Phi " << LogitTransformBack(mLogParam[0]) << " new Phi " << LogitTransformBack(mNewLogParam[0]) << endl;
//    cout<< "old Sigma2 " << exp(mLogParam[1]) << " new Sigma2 " << exp(mNewLogParam[1])<< endl;
//    cout<< "old " << dOldLogPosterior[0]<< " new " << dNewLogPosterior[0]<< endl;
//    cout<< "dLogU " << dLogU<<" dAccept " << dAccept<<  endl;

	if(dLogU<=dAccept)
	{
		sParam.dPhi[0]=LogitTransformBack(mNewLogParam[0]);
		sParam.dSigma2[0]=exp(mNewLogParam[1]);
//		 cout<< "XXXXXXXX Accepted XXXXXXXXX"<< endl;
	}
//	cout<< "XXXXXXXXXXXXXXXXX"<< endl;

//	cout<<"final Phi " << sParam.dPhi[0]<< endl;
//	cout<<"final Sigma2 " << sParam.dSigma2[0]<< endl;

	mLogParam[0]=LogitTransform(sParam.dPhi[0]);
	mLogParam[1]=log(sParam.dSigma2[0]);
	double dNumOfIter=(double) iNumOfIter;
	/* Covariance */
	if(iNumOfIter<=1)
	{
//		cout<<"first two iter END"<<endl;
//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;
		MatrixTrans( mLogParam, 2, 1, mLogParamTrans);
		MatrixMulti(mLogParam,mLogParam,2, 1,1,2, mParamParam);
		MatrixAdd(mCovar, mParamParam, 2,2, mCovar);
		MatrixAdd(mSum, mLogParam, 2,1, mSum);
//		cout<<"mSum[0] "<<mSum[0]<<endl;
//		cout<<"mSum[1] "<<mSum[1]<<endl;
//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

	}
	else if(iNumOfIter==2)
	{
//		cout<<"third iter END"<<endl;
		MatrixScalarMulti(mSigma,(dNumOfIter-1)/dNumOfIter,2,2,mCovar);

//		cout<<"mSigma[0] "<<mSigma[0]<<endl;
//		cout<<"mSigma[1] "<<mSigma[1]<<endl;
//		cout<<"mSigma[2] "<<mSigma[2]<<endl;
//		cout<<"mSigma[3] "<<mSigma[3]<<endl;
//
//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
//		cout<<"SumSum[0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum[1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum[2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum[3] "<<mSumSum[3]<<endl;
		MatrixScalarMulti(mSumSum,1/(dNumOfIter* dNumOfIter),2, 2, mSumSum);
//		cout<<"SumSum divided [0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum divided [1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum divided [2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum divided [3] "<<mSumSum[3]<<endl;
		MatrixAdd(mCovar, mSumSum, 2,2, mCovar);

//		cout<<"mSum[0] "<<mSum[0]<<endl;
//		cout<<"mSum[1] "<<mSum[1]<<endl;
//		cout<<"mCovar+SumSum[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum[3] "<<mCovar[3]<<endl;

		MatrixTrans( mLogParam, 2, 1, mLogParamTrans);
		MatrixMulti(mLogParam,mLogParamTrans,2, 1,1,2, mParamParam);
		MatrixScalarMulti(mParamParam,1/dNumOfIter,2, 2, mParamParam);
		MatrixAdd(mCovar, mParamParam, 2,2, mCovar);

//		cout<<"mParamParam[0] "<<mParamParam[0]<<endl;
//		cout<<"mParamParam[1] "<<mParamParam[1]<<endl;
//		cout<<"mParamParam[2] "<<mParamParam[2]<<endl;
//		cout<<"mParamParam[3] "<<mParamParam[3]<<endl;
//
//		cout<<"mCovar+SumSum+ParamParam[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[3] "<<mCovar[3]<<endl;
//
//
//		cout<<"mSum[0] "<<mSum[0]<<endl;
//		cout<<"mSum[1] "<<mSum[1]<<endl;
//		cout<<"mLogParam[0] "<<mLogParam[0]<<endl;
//		cout<<"mLogParam[1] "<<mLogParam[1]<<endl;
		mSum[0]=mSum[0]+ mLogParam[0];
		mSum[1]=mSum[1]+ mLogParam[1];

//		cout<<"mSum[0]+log"<<mSum[0]<<endl;
//		cout<<"mSum[1]+log"<<mSum[1]<<endl;

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
//		cout<<"SumSum[0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum[1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum[2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum[3] "<<mSumSum[3]<<endl;
		MatrixScalarMulti(mSumSum,1/(dNumOfIter*(dNumOfIter+1)),2, 2, mSumSum);
		MatrixSub(mCovar, mSumSum, 2,2, mCovar);
//		cout<<"SumSum divided [0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum divided [1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum divided [2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum divided [3] "<<mSumSum[3]<<endl;

//		cout<<"mCovar+SumSum+ParamParam+last[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[3] "<<mCovar[3]<<endl;


	}
	else
	{
//		cout<<"later iter END"<<endl;

//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

		MatrixScalarMulti(mCovar,(dNumOfIter-1)/dNumOfIter,2,2,mCovar);

//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
		MatrixScalarMulti(mSumSum,1/(dNumOfIter* dNumOfIter),2, 2, mSumSum);
		MatrixAdd(mCovar, mSumSum, 2,2, mCovar);

//		cout<<"mCovar+SumSum[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum[3] "<<mCovar[3]<<endl;

		MatrixTrans( mLogParam, 2, 1, mLogParamTrans);
//		cout<<"mLogParam[0] "<<mLogParam[0]<<endl;
//		cout<<"mLogParam[1] "<<mLogParam[1]<<endl;

		MatrixMulti(mLogParam,mLogParamTrans,2, 1,1,2, mParamParam);

//		cout<<"mParamParam[0] "<<mParamParam[0]<<endl;
//		cout<<"mParamParam[1] "<<mParamParam[1]<<endl;
//		cout<<"mParamParam[2] "<<mParamParam[2]<<endl;
//		cout<<"mParamParam[3] "<<mParamParam[3]<<endl;

		MatrixScalarMulti(mParamParam,1/( dNumOfIter),2, 2, mParamParam);
		MatrixAdd(mCovar, mParamParam, 2,2, mCovar);

//		cout<<"mCovar+SumSum+ParamParam[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[3] "<<mCovar[3]<<endl;

		mSum[0]=mSum[0]+LogitTransform(sParam.dPhi[0]);
		mSum[1]=mSum[1]+log(sParam.dSigma2[0]);

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
		MatrixScalarMulti(mSumSum,1/( dNumOfIter* (dNumOfIter+1)),2, 2, mSumSum);
		MatrixSub(mCovar, mSumSum, 2,2, mCovar);
//		cout<<"mCovar+SumSum+ParamParam+last[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[3] "<<mCovar[3]<<endl;


	}

	delete [] mChol;
	delete [] mSigma;
	delete [] mSumSum;
	delete [] mParamParam;
	delete [] mSumTrans;
	delete [] mLogParam;
	delete [] mLogParamTrans;
	delete [] mNewLogParam;
	delete [] mRandom;
	delete dOldLogPosterior;
	delete dNewLogPosterior;

}

void DrawPhi(struct AllParam sParam, int iNumOfObs, const gsl_rng * gsl_random_num)
{
	double* dSumSZt = new double;
	double * dSumCrossZt=new double;

	dSumSZt[0]=0;
	dSumCrossZt[0]=0;

	double dSigma20=0.1;
	double dPhi0=0.9;

	for(int i=1; i<iNumOfObs; i++)
	{
		dSumSZt[0]=dSumSZt[0]+ sParam.vX[i]*sParam.vX[i];
		dSumCrossZt[0]=dSumCrossZt[0]+ sParam.vX[i-1]*sParam.vX[i];
	}

	double dBetaOLS=dSumCrossZt[0]/dSumSZt[0];
	double dSigma2Star=1/(dSumSZt[0]/sParam.dSigma2[0] +1/dSigma20);
	double dBetaStar=dSigma2Star*(dSumSZt[0]*dBetaOLS/sParam.dSigma2[0]+ dPhi0/dSigma20);

	sParam.dPhi[0]=dBetaStar+gsl_ran_gaussian(gsl_random_num,sqrt(dSigma2Star));

	delete dSumSZt;
	delete dSumCrossZt;
}

void DrawSigma2(struct AllParam sParam, int iNumOfObs, const gsl_rng * gsl_random_num)
{
	double * dResS=new double;
	dResS[0]=0;
	for(int i=1; i<iNumOfObs; i++)
	{
		dResS[0]=dResS[0]+(sParam.vX[i]-sParam.vX[i-1]*sParam.dPhi[0])*(sParam.vX[i]-sParam.vX[i-1]*sParam.dPhi[0]);
	}
	int iNu=100;
	double dL=0.00001;

	sParam.dSigma2[0]=(iNu*dL+dResS[0])/gsl_ran_chisq (gsl_random_num, iNu+(iNumOfObs-1));

	delete dResS;
}

void DrawPhiSigmaFullConditional(struct AllParam sParam,  int iNumOfObs, const gsl_rng * gsl_random_num )
{

	 DrawPhi(sParam, iNumOfObs,gsl_random_num);
	 DrawSigma2(sParam, iNumOfObs,gsl_random_num);
}

void DrawN( struct AllParam sParam, double * vData,  int iNumOfObs , const gsl_rng * gsl_random_num, double *vZeroSkellamDens, double * vLogBessel )
{
	double* dCumU = new double;
	double* dLambda = new double;
	int* iTempN = new int;
	int* iNextBinom=new int;
	double* dU=new double;
	
	int MaxiTempN;
	MaxiTempN = 1000;
//	cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;


	for(int i=0; i<iNumOfObs; i++)
	{
//		cout<< "Draw N "<< i <<endl;
		dCumU[0]=0;
		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);

		iTempN[0]=fabs(vData[i]);

//		cout<<"XXXXXXXXXXXXXXXxx "<< i <<" XXXXXXXXXXXXXx"<<endl;
//		cout << "CumU nulla " << dCumU[0]<< endl;
//		cout << " sParam.dGamma[0] "<<  sParam.dGamma[0] << endl;
//		cout<< "sParam.vX[i] "<< sParam.vX[i]<<endl;
//		cout<< "sParam.vS[i] "<< sParam.vS[i]<<endl;
//		cout << " dLambda[0] "<<  dLambda[0] << endl;
//		cout << " vData[i] "<<  vData[i]<< endl;
//		cout << " vZeroSkellamDens "<<  vZeroSkellamDens[i]<< endl;

		if(vData[i]==0)
		{

//			dCumU[0]= (sParam.dGamma[0]*pow(2*dLambda[0], iTempN[0])*exp(-2*dLambda[0])/tgamma(iTempN[0]+1)+
//				      (1-sParam.dGamma[0])*pow(dLambda[0], iTempN[0])*exp(-2*dLambda[0])/(tgamma(0.5*(iTempN[0]+vData[i])+1)*
//				    		  tgamma(0.5*(iTempN[0]-vData[i])+1)))/vZeroSkellamDens[i];
			dCumU[0]=exp(log(sParam.dGamma[0])+iTempN[0]*log(2*dLambda[0])-2*dLambda[0]-lgamma(iTempN[0]+1)-log(vZeroSkellamDens[i]))+
					exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dLambda[0])-2*dLambda[0]-lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-log(vZeroSkellamDens[i]));
		}
		else
		{
//			dCumU[0]=((1-sParam.dGamma[0])*pow(dLambda[0], iTempN[0])*exp(-2*dLambda[0])/(tgamma(0.5*(iTempN[0]+vData[i])+1)*tgamma(0.5*(iTempN[0]-vData[i])+1)))/vZeroSkellamDens[i];
//			dCumU[0]=pow(dLambda[0], iTempN[0])/(tgamma(0.5*(iTempN[0]+vData[i])+1)*tgamma(0.5*(iTempN[0]-vData[i])+1)*exp(vLogBessel[i]));
			dCumU[0]=exp(iTempN[0]*log(dLambda[0])-lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogBessel[i]);
		}


		iNextBinom[0]=iTempN[0]+2;
		dU[0]=gsl_rng_uniform(gsl_random_num);

//		cout << "CumU " << dCumU[0] << " dU "<< dU[0]<< endl;

		// while(dCumU[0]<dU[0])
		while ((dCumU[0]<dU[0]) & (iTempN[0] < MaxiTempN))			
			
		{
			/*New N */
			iTempN[0]=iTempN[0]+1;

			if(vData[i]==0)
			{
				if(iNextBinom[0]==iTempN[0] )
				{
//					dCumU[0]=dCumU[0]+(sParam.dGamma[0]*pow(2*dLambda[0], iTempN[0])*exp(-2*dLambda[0])/tgamma(iTempN[0]+1)+
//					   (1-sParam.dGamma[0])*pow(dLambda[0], iTempN[0])*exp(-2*dLambda[0])/(tgamma(0.5*(iTempN[0]+vData[i])+1)*
//							   tgamma(0.5*(iTempN[0]-vData[i])+1)))/vZeroSkellamDens[i];

					dCumU[0]=dCumU[0]+exp(log(sParam.dGamma[0])+iTempN[0]*log(2*dLambda[0])-2*dLambda[0]-lgamma(iTempN[0]+1)-log(vZeroSkellamDens[i]))+
					exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dLambda[0])-2*dLambda[0]-lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-log(vZeroSkellamDens[i]));

					iNextBinom[0]=iTempN[0]+2;
				}
				else
				{
					//dCumU[0]=dCumU[0]+(sParam.dGamma[0]*pow(2*dLambda[0], iTempN[0])*exp(-2*dLambda[0])/tgamma(iTempN[0]+1))/vZeroSkellamDens[i];
					//dCumU[0]=dCumU[0]+sParam.dGamma[0]*exp( iTempN[0]*log(2*dLambda[0])-lgamma(iTempN[0]+1)- 2*dLambda[0])/vZeroSkellamDens[i];
					dCumU[0]=dCumU[0]+exp(log(sParam.dGamma[0])+iTempN[0]*log(2*dLambda[0])-2*dLambda[0]-lgamma(iTempN[0]+1)-log(vZeroSkellamDens[i]));
				}

			}
			else
			{
				if(iNextBinom[0]==iTempN[0] )
				{
//					dCumU[0]=dCumU[0]+((1-sParam.dGamma[0])*pow(dLambda[0], iTempN[0])*exp(-2*dLambda[0])/(tgamma(0.5*(iTempN[0]+vData[i])+1)*
//							tgamma(0.5*(iTempN[0]-vData[i])+1)))/vZeroSkellamDens[i];
//					dCumU[0]=dCumU[0]+pow(dLambda[0], iTempN[0])/(tgamma(0.5*(iTempN[0]+vData[i])+1)*tgamma(0.5*(iTempN[0]-vData[i])+1)*exp(vLogBessel[i]));
					dCumU[0]=dCumU[0]+exp(iTempN[0]*log(dLambda[0])-lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogBessel[i]);

//					cout << "iTempN[0]*log(dLambda[0]) " << iTempN[0]*log(dLambda[0])<< " log(vBessel[i]) "<< vLogBessel[i]<< endl;
					iNextBinom[0]=iTempN[0]+2;
				}
			}
//			if(iTempN[0]<5000)
//			{
//			cout << "vData[i] " << vData[i]<< " i " << i << endl;
//			cout << "CumU " << dCumU[0] << " dU "<< dU[0]<< endl;
//			cout << "Contirbuiton " <<exp(iTempN[0]*log(dLambda[0])-lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogBessel[i])<<endl;
//			cout << "iNextBinom[0]  " << iNextBinom[0]<< " dLambda[0] "<< dLambda[0]<< endl;
//
//			}
		}

		sParam.vN[i]=iTempN[0];



		if(!isthisfinite( sParam.vN[i]))
		{
			cout <<"N enter"<<endl;
			double * vNError = new double [11];

			vNError[0]=i;
			vNError[1]=sParam.vN[i];
			vNError[2]=dCumU[0];
			vNError[3]=dU[0];
			vNError[4]=dLambda[0];
			vNError[5]=vLogBessel[i];
			vNError[6]=vZeroSkellamDens[i];
			vNError[7]=vData[i];
			vNError[8]=sParam.dMu[0];
			vNError[9]=sParam.vS[i];
			vNError[10]=sParam.vX[i];
			cout <<"N write"<<endl;
			string sFile="vNError.csv" ;
			WriteOutDoubleArray( vNError , 11, 1, sFile);
			cout <<"N del"<<endl;
			delete [] vNError;
			exit(1);

		}



		// cout<< "sParam.vN[i] "<<sParam.vN[i]<<endl;
//		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX "<< endl;
	}



	delete dCumU;
	delete dLambda;
	delete iTempN;
	delete iNextBinom;
	delete dU;

}

void DrawTau(struct AllParam sParam,int iNumOfObs , const gsl_rng * gsl_random_num)
{


	for(int i=0; i<iNumOfObs; i++)
	{
//		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX "<< endl;
//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
//		cout << "sParam.vN[i] " << sParam.vN[i]<< endl;
//		cout<< "sParam.vX[i] "<< sParam.vX[i]<<endl;
		double dU2=gsl_rng_uniform(gsl_random_num);
		while(dU2==0)
		{
			dU2=gsl_rng_uniform(gsl_random_num);
		}

		/* Drawing Beta(N,1)*/
		if(sParam.vN[i]>0)
		{

			sParam.vTau2[i]=pow(dU2, (double ) 1.0/sParam.vN[i]);




//			cout<< " 1/sParam.vN[i] "<< (double) 1/sParam.vN[i] << endl;
//			cout << "dU2 " << dU2<< endl;
//			cout << "sParam.vTau2[i] " << sParam.vTau2[i]<< endl;
		}
		else
		{
			sParam.vTau2[i]=0;
		}
		/* Drawing exponential with intensity 2 lambda*/

		double dU1=gsl_rng_uniform(gsl_random_num);
		sParam.vTau1[i]=1-log(1-dU1)/(2*exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]))-sParam.vTau2[i];


		if(!isthisfinite(-log(sParam.vTau2[i])) & (sParam.vN[i]>0))
		{


			string sFile="vTau2Error.csv";
			WriteOutDoubleArray( sParam.vTau2, iNumOfObs, 1, sFile);
			sFile="dITau2Error.csv";
			WriteOutIntArray( &i, 1, 1, sFile);
			sFile="dUTau2Error.csv";
			WriteOutDoubleArray( &dU2, 1, 1, sFile);
			sFile="vNErrorTau.csv";
			WriteOutIntArray( sParam.vN, iNumOfObs, 1, sFile);
			sFile="vTau1Error.csv";
			WriteOutDoubleArray( sParam.vTau1, iNumOfObs, 1, sFile);
			exit(1);
		}

//		cout << "dU1 " << dU1<< endl;
//
//		cout << "sParam.vTau1[i] " << sParam.vTau1[i]<< endl;


	}
//	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX "<< endl;
}

double digammaOnHost(int n)
{
	/* Digamma function for integers using the recursion psi(n)=psi(n-1)+1/(n-1) */
	double dTemp=-0.57721566490153286060651209008240243104215933593992359880576;
	for(int i=2; i<=n;i++)
	{

		dTemp=dTemp+(double) 1/(i-1);
	}
	return dTemp;

}

double trigammaOnHost(int n)
{
	/* Trigamma function for integers using the recursion psi_1(n)=psi_1(n-1)-1/(n-1)² */
	double dTemp=1.6449340668482264364724151666460251892189499012067984;
	for(int i=2; i<=n;i++)
	{

		dTemp=dTemp-(double) 1/((i-1)*(i-1));
	}
	return dTemp;
}

void DrawAuxYandAuxH(struct AllParam sParam, int iNumOfObs , const gsl_rng * gsl_random_num)
{
	int iNumOfComp;
	double* vWeights =new double [10];
	double* vMeans = new double [10];
	double* vVariances= new double [10];
	double* dU=new double;
	double* dCumP=new double[10];
	double* dLogLambda=new double;

	for(int k=0; k<iNumOfObs; k++)
	{

//		cout<<"XXXXXXXXXXXXXXXX  Start XXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;
//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
//		cout << "sParam.vN[i] " << sParam.vN[k]<< endl;
//		cout<< "sParam.vX[i] "<< sParam.vX[k]<<endl;
//		cout << "sParam.vTau1[k] " << sParam.vTau1[k]<<endl;
//		cout << "sParam.vTau2[k] " << sParam.vTau2[k]<<endl;

		dLogLambda[0]=sParam.dMu[0]+sParam.vS[k]+sParam.vX[k];

		/* Mixture approximation for tau1 */
		double mu = - digammaOnHost(1);
		double sigma = sqrt(trigammaOnHost(1));

		/* Weights */
		vWeights[0]=0.00396984425; vWeights[1]=0.0396244597; vWeights[2]=0.16776747; vWeights[3]=0.147036501; vWeights[4]=0.125306271;
		vWeights[5]=0.1014852; vWeights[6]=0.103758531; vWeights[7]=0.115972617; vWeights[8]=0.107065659; vWeights[9]=0.0880134482;

		/* Means */
		vMeans[0]=3.55887454*sigma+mu; vMeans[1]=2.11415904*sigma+mu; vMeans[2]= 0.968124631*sigma+mu; vMeans[3]= 0.51537638*sigma+mu; vMeans[4]= 0.145465449*sigma+mu;
		vMeans[5]=-0.145346445*sigma+mu; vMeans[6]= -0.416660312*sigma+mu; vMeans[7]=-0.689002855*sigma+mu; vMeans[8]= -0.974965634*sigma+mu; vMeans[9]= -1.27310004*sigma+mu;

		/* Variances */
		vVariances[0]=2.62603032*sigma*sigma; vVariances[1]=1.21263644*sigma*sigma; vVariances[2]=0.66586521*sigma*sigma; vVariances[3]=0.256650604*sigma*sigma; vVariances[4]= 0.120071142*sigma*sigma;
		vVariances[5]=0.0649909219*sigma*sigma; vVariances[6]=0.0473513798*sigma*sigma; vVariances[7]= 0.046639443*sigma*sigma; vVariances[8]= 0.0576541602*sigma*sigma; vVariances[9]= 0.0888536903*sigma*sigma;

//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//
//		}
//
//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//		}
//
//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//		}

		/* Draw indicator form the mixture approximation for tau1*/
//		dCumP[0]=vWeights[0]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0])/sqrt(vVariances[0]);
		dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0]);
//		dCumP[0]=exp(log(vWeights[0])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0], sqrt(vVariances[0]))));
		for(int i=1; i<10;i++)
		{
//			dCumP[i]=dCumP[i-1]+vWeights[i]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i])/sqrt(vVariances[i]);
			dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i]);
//			dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i], sqrt(vVariances[i]))));
//			cout<<"XX  contribution "<< vWeights[i]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i])/sqrt(vVariances[i])<<endl;
//			cout<<"XX  dCumP[i-1] "<< dCumP[i-1]<<endl;
//			cout<<"XX  dCumP[i] "<< dCumP[i]<<endl;

		}

//		for(int l=0; l<10; l++)
//		{
//			cout<<"1 dCumP[i]/dCumP[iNumOfComp-1] at "<<l << " is " << dCumP[l]/dCumP[9]<< endl;
//
//		}


		int i=0;
		dU[0]=gsl_rng_uniform(gsl_random_num);
//		cout<<"dU 1 " << dU[0]<<endl;
		while((dCumP[i]/dCumP[10-1])<dU[0])
		{
			i=i+1;
//			cout<<"i "<< i << endl;
		}

		/*Saving AuxY and AuxH for Tau1  */

//		cout<<"-log(sParam.vTau1[k]) "<<-log(sParam.vTau1[k])<<endl;
//		cout <<" -log 2 "<<-log(2)<<endl;
//		cout<<"-vMeans[i] "<<-vMeans[i] <<endl;
		sParam.mAuxY[2*k]=-log(sParam.vTau1[k])-log(2)-vMeans[i];//-sParam.vS[k] ;
		sParam.mAuxH[4*k]=vVariances[i];


		if(sParam.vN[k]>0)
		{
			/* Mixture approximation for tau2 */
			MixtureLogGamma( sParam.vN[k], &iNumOfComp,  vWeights , vMeans ,vVariances);
			/* Draw indicator form the mixture approximation for tau2*/
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//
//			}
//
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//			}
//
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//			}

//				dCumP[0]=vWeights[0]*exp(-0.5*pow((-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0]),2) /vVariances[0])/sqrt(vVariances[0]);
				dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0]);
//				dCumP[0]=exp(log(vWeights[0])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0], sqrt(vVariances[0]))));
				for(int i=1; i<iNumOfComp;i++)
				{
//					dCumP[i]=dCumP[i-1]+vWeights[i]*exp(-0.5*pow((-log(sParam.vTau2[k])-log( 2)-dLogLambda[0]-vMeans[i]),2) /vVariances[i])/sqrt(vVariances[i]);
					dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i]);
//					dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[i], sqrt(vVariances[i]))));

//					cout<<"XX  contribution "<< vWeights[i]*exp(-0.5*pow((-log(sParam.vTau2[k])-log( 2)-dLogLambda[0]-vMeans[i]),2) /vVariances[i])/sqrt(vVariances[i])<<endl;
//					cout<<"XX  dCumP[i-1] "<< dCumP[i-1]<<endl;
//					cout<<"XX  dCumP[i] "<< dCumP[i]<<endl;
				}


//				for(int l=0; l<iNumOfComp; l++)
//				{
//					cout<<"sParam.vN[k] "<<sParam.vN[k]<<" iNumOfComp "<< iNumOfComp << endl;
//					cout<<"2 dCumP[i]/dCumP[iNumOfComp-1] at "<<l << " is " << dCumP[l]/dCumP[iNumOfComp-1]<< endl;
//
//				}

				int j=0;
				dU[0]=gsl_rng_uniform(gsl_random_num);
//				cout<<"dU 2 " << dU[0]<<endl;
				while(dCumP[j]/dCumP[iNumOfComp-1]<dU[0])
				{
					j=j+1;
//					cout<<"j "<< j << endl;
				}

				/*Saving AuxY and AuxH for Tau2 */


				if(!isthisfinite( -log(sParam.vTau2[k])-log(2)-vMeans[j]))
										{
											double * vYErrorAux = new double [4];

											vYErrorAux[0]=k;
											vYErrorAux[1]=sParam.vTau2[k];
											vYErrorAux[1]=-log(sParam.vTau2[k]);
											vYErrorAux[2]=vMeans[j];



											string sFile="vYErrorAux.csv";
											WriteOutDoubleArray( vYErrorAux , 4, 1, sFile);

											delete [] vYErrorAux;
											cout<<"AuxY Error"<<endl;
											exit(1);

										}


				if(!isthisfinite( 1.0/vVariances[j]))
								{
									double * vVarErrorAux = new double [iNumOfComp+2];

									for(int t=0; t<iNumOfComp;t++)
									{
										vVarErrorAux[t]=vVariances[t];
									}

									vVarErrorAux[iNumOfComp]=sParam.vN[k];
									vVarErrorAux[iNumOfComp+1]=k;
									string sFile="vVarErrorAux.csv";
									WriteOutDoubleArray( vVarErrorAux , iNumOfComp+2, 1, sFile);

									cout<<"AuxH Error"<<endl;
									delete [] vVarErrorAux;
									exit(1);
					}


				sParam.mAuxY[2*k+1]= -log(sParam.vTau2[k])-log(2)-vMeans[j];//-sParam.vS[k];
				sParam.mAuxH[4*k+3]=vVariances[j] ;

		}
		else
		{
			sParam.mAuxY[2*k+1]= 0;
			sParam.mAuxH[4*k+3]=0;

		}


//		cout<<"AuxY 0 at "<<k << " is " <<sParam.mAuxY[2*k] << endl;
//		cout<<"AuxY 1 at "<<k << " is " <<sParam.mAuxY[2*k+1] << endl;
//		cout<<"AuxH 0 at "<<k << " is " <<sParam.mAuxH[4*k] << endl;
//		cout<<"AuxH 1 at "<<k << " is " <<sParam.mAuxH[4*k+3] << endl;
//
//		cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;

		sParam.mAuxH[4*k+1]=0;
		sParam.mAuxH[4*k+2]=0;


	}

	delete [] vWeights;
	delete [] vMeans;
	delete [] vVariances;
	delete [] dCumP;

	delete dU;
	delete dLogLambda;

}

void SimulateAuxYandAuxH(struct AllParam sParam, int iNumOfObs , const gsl_rng * gsl_random_num)
{
	int iNumOfComp;
	double* vWeights =new double [10];
	double* vMeans = new double [10];
	double* vVariances= new double [10];
	double* dU=new double;
	double* dCumP=new double[10];
	double* dLogLambda=new double;

	for(int k=0; k<iNumOfObs; k++)
	{

//		cout<<"XXXXXXXXXXXXXXXX  Start XXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;
//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
//		cout << "sParam.vN[i] " << sParam.vN[k]<< endl;
//		cout<< "sParam.vX[i] "<< sParam.vX[k]<<endl;
//		cout << "sParam.vTau1[k] " << sParam.vTau1[k]<<endl;
//		cout << "sParam.vTau2[k] " << sParam.vTau2[k]<<endl;

		dLogLambda[0]=sParam.dMu[0]+sParam.vS[k]+sParam.vX[k];

		/* Mixture approximation for tau1 */
		double mu = - digammaOnHost(1);
		double sigma = sqrt(trigammaOnHost(1));

		/* Weights */
		vWeights[0]=0.00396984425; vWeights[1]=0.0396244597; vWeights[2]=0.16776747; vWeights[3]=0.147036501; vWeights[4]=0.125306271;
		vWeights[5]=0.1014852; vWeights[6]=0.103758531; vWeights[7]=0.115972617; vWeights[8]=0.107065659; vWeights[9]=0.0880134482;

		/* Means */
		vMeans[0]=3.55887454*sigma+mu; vMeans[1]=2.11415904*sigma+mu; vMeans[2]= 0.968124631*sigma+mu; vMeans[3]= 0.51537638*sigma+mu; vMeans[4]= 0.145465449*sigma+mu;
		vMeans[5]=-0.145346445*sigma+mu; vMeans[6]= -0.416660312*sigma+mu; vMeans[7]=-0.689002855*sigma+mu; vMeans[8]= -0.974965634*sigma+mu; vMeans[9]= -1.27310004*sigma+mu;

		/* Variances */
		vVariances[0]=2.62603032*sigma*sigma; vVariances[1]=1.21263644*sigma*sigma; vVariances[2]=0.66586521*sigma*sigma; vVariances[3]=0.256650604*sigma*sigma; vVariances[4]= 0.120071142*sigma*sigma;
		vVariances[5]=0.0649909219*sigma*sigma; vVariances[6]=0.0473513798*sigma*sigma; vVariances[7]= 0.046639443*sigma*sigma; vVariances[8]= 0.0576541602*sigma*sigma; vVariances[9]= 0.0888536903*sigma*sigma;

//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//
//		}
//
//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//		}
//
//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//		}

		/* Draw indicator form the mixture approximation for tau1*/
//		dCumP[0]=vWeights[0]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0])/sqrt(vVariances[0]);
		dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0]);
//		dCumP[0]=exp(log(vWeights[0])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0], sqrt(vVariances[0]))));
		for(int i=1; i<10;i++)
		{
//			dCumP[i]=dCumP[i-1]+vWeights[i]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i])/sqrt(vVariances[i]);
			dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i]);
//			dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i], sqrt(vVariances[i]))));
//			cout<<"XX  contribution "<< vWeights[i]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i])/sqrt(vVariances[i])<<endl;
//			cout<<"XX  dCumP[i-1] "<< dCumP[i-1]<<endl;
//			cout<<"XX  dCumP[i] "<< dCumP[i]<<endl;

		}

//		for(int l=0; l<10; l++)
//		{
//			cout<<"1 dCumP[i]/dCumP[iNumOfComp-1] at "<<l << " is " << dCumP[l]/dCumP[9]<< endl;
//
//		}


		int i=0;
		dU[0]=gsl_rng_uniform(gsl_random_num);
//		cout<<"dU 1 " << dU[0]<<endl;
		while((dCumP[i]/dCumP[10-1])<dU[0])
		{
			i=i+1;
//			cout<<"i "<< i << endl;
		}

		/*Saving AuxY and AuxH for Tau1  */

//		cout<<"-log(sParam.vTau1[k]) "<<-log(sParam.vTau1[k])<<endl;
//		cout <<" -log 2 "<<-log(2)<<endl;
//		cout<<"-vMeans[i] "<<-vMeans[i] <<endl;
		sParam.mAuxY[2*k]=sParam.dMu[0]+sParam.vS[k]+sParam.vX[k]+gsl_ran_gaussian(gsl_random_num,pow(vVariances[i],0.5));//-sParam.vS[k] ;
		sParam.mAuxH[4*k]=vVariances[i];


		if(sParam.vN[k]>0)
		{
			/* Mixture approximation for tau2 */
			MixtureLogGamma( sParam.vN[k], &iNumOfComp,  vWeights , vMeans ,vVariances);
			/* Draw indicator form the mixture approximation for tau2*/
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//
//			}
//
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//			}
//
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//			}

//				dCumP[0]=vWeights[0]*exp(-0.5*pow((-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0]),2) /vVariances[0])/sqrt(vVariances[0]);
				dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0]);
//				dCumP[0]=exp(log(vWeights[0])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0], sqrt(vVariances[0]))));
				for(int i=1; i<iNumOfComp;i++)
				{
//					dCumP[i]=dCumP[i-1]+vWeights[i]*exp(-0.5*pow((-log(sParam.vTau2[k])-log( 2)-dLogLambda[0]-vMeans[i]),2) /vVariances[i])/sqrt(vVariances[i]);
					dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i]);
//					dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[i], sqrt(vVariances[i]))));

//					cout<<"XX  contribution "<< vWeights[i]*exp(-0.5*pow((-log(sParam.vTau2[k])-log( 2)-dLogLambda[0]-vMeans[i]),2) /vVariances[i])/sqrt(vVariances[i])<<endl;
//					cout<<"XX  dCumP[i-1] "<< dCumP[i-1]<<endl;
//					cout<<"XX  dCumP[i] "<< dCumP[i]<<endl;
				}


//				for(int l=0; l<iNumOfComp; l++)
//				{
//					cout<<"sParam.vN[k] "<<sParam.vN[k]<<" iNumOfComp "<< iNumOfComp << endl;
//					cout<<"2 dCumP[i]/dCumP[iNumOfComp-1] at "<<l << " is " << dCumP[l]/dCumP[iNumOfComp-1]<< endl;
//
//				}

				int j=0;
				dU[0]=gsl_rng_uniform(gsl_random_num);
//				cout<<"dU 2 " << dU[0]<<endl;
				while(dCumP[j]/dCumP[iNumOfComp-1]<dU[0])
				{
					j=j+1;
//					cout<<"j "<< j << endl;
				}

				/*Saving AuxY and AuxH for Tau2 */

				//sParam.mAuxY[2*k+1]= log(2)+sParam.dMu[0]+sParam.vS[k]+sParam.vX[k]+vMeans[j]+gsl_ran_gaussian(gsl_random_num,pow(vVariances[j],0.5));//-sParam.vS[k];
				sParam.mAuxY[2*k+1]= sParam.dMu[0]+sParam.vS[k]+sParam.vX[k]+gsl_ran_gaussian(gsl_random_num,pow(vVariances[j],0.5));//-sParam.vS[k];
				sParam.mAuxH[4*k+3]=vVariances[j] ;

		}
		else
		{
			sParam.mAuxY[2*k+1]= 0;
			sParam.mAuxH[4*k+3]=0;

		}


//		cout<<"AuxY 0 at "<<k << " is " <<sParam.mAuxY[2*k] << endl;
//		cout<<"AuxY 1 at "<<k << " is " <<sParam.mAuxY[2*k+1] << endl;
//		cout<<"AuxH 0 at "<<k << " is " <<sParam.mAuxH[4*k] << endl;
//		cout<<"AuxH 1 at "<<k << " is " <<sParam.mAuxH[4*k+3] << endl;
//
//		cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;

		sParam.mAuxH[4*k+1]=0;
		sParam.mAuxH[4*k+2]=0;


	}

	delete [] vWeights;
	delete [] vMeans;
	delete [] vVariances;
	delete [] dCumP;

	delete dU;
	delete dLogLambda;

}


void CaculateDNB(const double * vData, int iNumOfObs,  struct AllParam sParam,double * vDNB, double * vZeroDNBDens,
		double * vIndicator)
{
	double* dLambda=new double;
	double* dLambdaRatio=new double;
	double* dNuRatio=new double;
	int * iAbsX=new int;

	for(int i=0; i<iNumOfObs; i++)
	{
//		cout <<"log lambda "<< <<endl;
		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);
		dLambdaRatio[0]=dLambda[0]/(dLambda[0]+sParam.dNu[0]);
		//dNuRatio[0]=sParam.dNu[0]/(dLambda[0]+sParam.dNu[0]);
		dNuRatio[0]=1-dLambdaRatio[0];
		iAbsX[0]=abs(vData[i]);

//		cout << "Hyper original Star "<<endl;
//		cout << "Hyper "<<gsl_sf_hyperg_2F1(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2))<<endl;
//		cout << "Hyper original End "<<endl;
//		cout << "Hyper Star "<<endl;
//		cout << "Hyper "<<gsl_sf_hyperg_2F1_renorm(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2))<<endl;
//		cout << "Hyper End "<<endl;

//		cout <<"log DNB "<<2*sParam.dNu[0]*log(dNuRatio[0])+iAbsX[0]*log(dLambdaRatio[0])+lgamma(sParam.dNu[0]+iAbsX[0])-lgamma(sParam.dNu[0])
//						+log(gsl_sf_hyperg_2F1_renorm(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2)))<<endl;

//		vDNB[i]=exp(2*sParam.dNu[0]*log(dNuRatio[0])+iAbsX[0]*log(dLambdaRatio[0])+lgamma(sParam.dNu[0]+iAbsX[0])-lgamma(sParam.dNu[0])
//				+log(gsl_sf_hyperg_2F1_renorm(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2))));

//		cout <<"vData[i] "<< vData[i] <<endl;
//		cout <<"Log Lambda "<<sParam.dMu[0]+sParam.vX[i]<<endl;
		gsl_sf_result result;
		int status = gsl_sf_hyperg_2F1_e(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2) ,&result);
		if (status) { /* an error occurred */
		/* status value specifies the type of error */
			double * vDNBError = new double [4];

			vDNBError[0]=sParam.dNu[0]+iAbsX[0];
			vDNBError[1]=sParam.dNu[0];
			vDNBError[2]=iAbsX[0]+1;
			vDNBError[3]=pow(dLambdaRatio[0],2);



			string sFile="vDNBError.csv";
			WriteOutDoubleArray( vDNBError , 4, 1, sFile);

			delete [] vDNBError;

			vDNB[i]=exp(2*sParam.dNu[0]*log(dNuRatio[0])+iAbsX[0]*log(dLambdaRatio[0])+lgamma(sParam.dNu[0]+iAbsX[0])-lgamma(sParam.dNu[0])-lgamma(iAbsX[0]+1)
											+log(HyperGeo(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2))));

		}
		else{

			vDNB[i]=exp(2*sParam.dNu[0]*log(dNuRatio[0])+iAbsX[0]*log(dLambdaRatio[0])+lgamma(sParam.dNu[0]+iAbsX[0])-lgamma(sParam.dNu[0])-lgamma(iAbsX[0]+1)
								+log(result.val));
		}

//		vDNB[i]=exp(2*sParam.dNu[0]*log(dNuRatio[0])+iAbsX[0]*log(dLambdaRatio[0])+lgamma(sParam.dNu[0]+iAbsX[0])-lgamma(sParam.dNu[0])-lgamma(iAbsX[0]+1)
//					+log(gsl_sf_hyperg_2F1(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2))));

		if(vData[i]==0)
		{
			vIndicator[i]=1;
		}
		else
		{
			vIndicator[i]=0;
		}

//		cout << "vIndicator "<<vIndicator[i]<<endl;
		vZeroDNBDens[i]=sParam.dGamma[0]*vIndicator[i]+(1-sParam.dGamma[0])*vDNB[i] ;
	}

	delete dLambda;
	delete dLambdaRatio;
	delete dNuRatio;
	delete iAbsX;
}

void DrawZ(struct AllParam sParam, double* vData,  int iNumOfObs ,const gsl_rng * gsl_random_num, double* vIndicator, double *vLogSkellamZ)
{
	double* dZ1=new double;
	double* dZ2=new double;
	double* dLogSkellamNew=new double;
	double* dLogSkellamOld=new double;
	double* dLogAccept=new double;

	double* dLambda=new double;

	for(int i=0; i<iNumOfObs; i++)
	{

//		cout<< "log Lambda in Z at "<< i <<" is "<< sParam.dMu[0]+sParam.vS[i]+sParam.vX[i] <<endl;
		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);
//		cout << "log lambda utan"<<endl;
		/* Draw Z1 and Z2 */
		dZ1[0]=gsl_ran_gamma (gsl_random_num, sParam.dNu[0], 1/sParam.dNu[0]);
		dZ2[0]=gsl_ran_gamma (gsl_random_num, sParam.dNu[0], 1/sParam.dNu[0]);
		/* Calculate acceptance probability */
		dLogSkellamNew[0]= -dLambda[0]*dZ1[0]-dLambda[0]*dZ2[0] +2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5 )+
				0.5*vData[i]*log(dZ1[0]/dZ2[0]);
//		cout<< "XXXXXXX "<< i << " XXXXXXXXXX"<<endl;
//		cout<< "Data[i] "<< vData[i]<<endl;
//		cout << "new Skellam start "<<endl;
		if(1e-50>2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5))
		{
//			cout<<"small"<<endl;
			dLogSkellamNew[0]=dLogSkellamNew[0]-2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5 )+abs(vData[i]) * log(dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5))-lgamma(abs(vData[i])+1);
		}
		else if(700<2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5))
		{
//			cout<<"large"<<endl;
			dLogSkellamNew[0]=dLogSkellamNew[0]-2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5 )+LogModifiedBesselFirstKindLargeX(abs(vData[i]), 2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5));
		}
		else
		{
//			cout << "middle "<<2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5)<<endl;
//			cout<<"2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5) "<< 2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5)<<endl;
//			cout<<"Bessel start"<<endl;
//			cout<<"Bessel "<< gsl_sf_bessel_In(abs(vData[i]), 100)<<endl;
//			cout<<"Bessel end"<<endl;
	//		dLogSkellamNew[0]=dLogSkellamNew[0]+exp(-2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5))*gsl_sf_bessel_In(abs(vData[i]), 2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5));
//			gsl_sf_result* result;
//			gsl_sf_bessel_In_scaled_e (abs(vData[i]), 2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5), result);
//			cout<< "error" <<result->err <<endl;
//			cout<< "val" <<result->val <<endl;
			dLogSkellamNew[0]=dLogSkellamNew[0]+log(gsl_sf_bessel_In_scaled(abs(vData[i]),2*dLambda[0]*pow(dZ1[0]*dZ2[0], 0.5)));

		}
//		cout << "new Skellam end "<<endl;

		dLogSkellamOld[0]= -dLambda[0]*sParam.vZ1[i]-dLambda[0]*sParam.vZ2[i] +2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5 )+
						0.5*vData[i]*log(sParam.vZ1[i]/sParam.vZ2[i]);

//		cout << "old Skellam start "<<endl;
		if(1e-50>2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5))
		{
			dLogSkellamOld[0]=dLogSkellamOld[0]-2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5 )+abs(vData[i]) * log(dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5))-lgamma(abs(vData[i])+1);
		}
		else if(700<2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5))
		{
			dLogSkellamOld[0]=dLogSkellamOld[0]-2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5 )+LogModifiedBesselFirstKindLargeX(abs(vData[i]), 2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5));
		}
		else
		{
			dLogSkellamOld[0]=dLogSkellamOld[0]+log(gsl_sf_bessel_In_scaled(abs(vData[i]),2*dLambda[0]*pow(sParam.vZ1[i]*sParam.vZ2[i], 0.5)));
		}
//		cout << "old Skellam end "<<endl;

		dLogAccept[0]=fmin(log(sParam.dGamma[0]*vIndicator[i]+(1-sParam.dGamma[0])*exp(dLogSkellamNew[0])) -log(sParam.dGamma[0]*vIndicator[i]+(1-sParam.dGamma[0])*exp(dLogSkellamOld[0])),0);



		/* Reject step  and save Skellam*/
		double dLogU=log(gsl_rng_uniform(gsl_random_num));

		if((dLogU<=dLogAccept[0])& (isthisfinite(dZ1[0]))&  (isthisfinite(dZ2[0])) )
		{
			sParam.vZ1[i]=dZ1[0];
			sParam.vZ2[i]=dZ2[0];
			vLogSkellamZ[i]=log(sParam.dGamma[0]*vIndicator[i]+(1-sParam.dGamma[0])*exp(dLogSkellamNew[0]));
		}
		else
		{
			vLogSkellamZ[i]=log(sParam.dGamma[0]*vIndicator[i]+(1-sParam.dGamma[0])*exp(dLogSkellamOld[0]));
		}

//		cout<< "vData[i] "<<vData[i]<<endl;
//		cout<< "dZ1[0] " << dZ1[0] <<endl;
//		cout<< "dZ2[0] " << dZ2[0] <<endl;
//		cout<< "dLogSkellamNew[0] " << dLogSkellamNew[0] <<endl;
//		cout<< "dLogSkellamOld[0] " << dLogSkellamOld[0] <<endl;
//		cout<< "dLogAccept[0] " << dLogAccept[0] <<endl;
//		cout<< "sParam.dGamma[0]*vIndicator[i] " << sParam.dGamma[0]*vIndicator[i] <<endl;
//		cout<< "vLogSkellamZ[i] " << vLogSkellamZ[i]<<endl;

	}

	delete dLambda;
	delete dZ1;
	delete dZ2;
	delete dLogSkellamNew;
	delete dLogSkellamOld;
	delete dLogAccept;

}

void DrawNDNB( struct AllParam sParam, double * vData,  int iNumOfObs , const gsl_rng * gsl_random_num, double *vLogSkellamZ )
{
	double* dCumU = new double;
	double* dLambda = new double;
	double* dZSum=new double;
	double* dZ1Ratio=new double;
	double* dZ2Ratio=new double;
	double* dZLambda=new double;
	int* iTempN = new int;
	 int* iNextBinom=new int;
	double* dU=new double;

	int MaxiTempN;
	MaxiTempN = 1000;
//	cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;

	for(int i=0; i<iNumOfObs; i++)
	{
//		cout<< "Draw N "<< i <<endl;
		dCumU[0]=0;
		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);
		iTempN[0]=fabs(vData[i]);
		dZSum[0]=sParam.vZ1[i]+sParam.vZ2[i];
		dZ1Ratio[0]=sParam.vZ1[i]/dZSum[0];
		dZ2Ratio[0]=sParam.vZ2[i]/dZSum[0];
		dZLambda[0]=dZSum[0]*dLambda[0];
//		cout<<"vData[i] "<< vData[i]<<endl;
//		cout<<"dZSum[0] "<<dZSum[0]<<endl;
//		cout<<"dZ1Ratio[0] "<<dZ1Ratio[0]<<endl;
//		cout<<"dZ2Ratio[0] "<<dZ2Ratio[0]<<endl;
//		cout<<"dZLambda[0] "<<dZLambda[0]<<endl;

		if(vData[i]==0)
		{
			dCumU[0]=exp(log(sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]-lgamma(iTempN[0]+1)-vLogSkellamZ[i])+
					exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]+0.5*(iTempN[0]+vData[i])*log(dZ1Ratio[0])+0.5*(iTempN[0]-vData[i])*log(dZ2Ratio[0])
					     -lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogSkellamZ[i]);
		}
		else
		{
			dCumU[0]=exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]+0.5*(iTempN[0]+vData[i])*log(dZ1Ratio[0])+0.5*(iTempN[0]-vData[i])*log(dZ2Ratio[0])
				     -lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogSkellamZ[i]);
		}


		iNextBinom[0]=iTempN[0]+2;
		dU[0]=gsl_rng_uniform(gsl_random_num);
		
		 
		// cout << "CumU " << dCumU[0] << " dU "<< dU[0]<< endl;
		// while(dCumU[0]<dU[0])
		while ((dCumU[0]<dU[0]) & (iTempN[0] < MaxiTempN))			
//		while(dCumU[0]<1)
		{
			/*New N */
			iTempN[0]=iTempN[0]+1;

			if(vData[i]==0)
			{
				if(iNextBinom[0]==iTempN[0] )
				{
					dCumU[0]=dCumU[0]+exp(log(sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]-lgamma(iTempN[0]+1)-vLogSkellamZ[i])+
							exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]+0.5*(iTempN[0]+vData[i])*log(dZ1Ratio[0])+0.5*(iTempN[0]-vData[i])*log(dZ2Ratio[0])
							     -lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogSkellamZ[i]);
					iNextBinom[0]=iTempN[0]+2;
				}
				else
				{

					dCumU[0]=dCumU[0]+exp(log(sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]-lgamma(iTempN[0]+1)-vLogSkellamZ[i]);
				}

			}
			else
			{
				if(iNextBinom[0]==iTempN[0] )
				{
					dCumU[0]=dCumU[0]+exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]+0.5*(iTempN[0]+vData[i])*log(dZ1Ratio[0])+0.5*(iTempN[0]-vData[i])*log(dZ2Ratio[0])
						     -lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogSkellamZ[i]);
					iNextBinom[0]=iTempN[0]+2;
				}
			}


//			if(iTempN[0]>1000)
//			{
//				cout << "contribution "<<exp(log(1-sParam.dGamma[0])+iTempN[0]*log(dZLambda[0])-dZLambda[0]+0.5*(iTempN[0]+vData[i])*log(dZ1Ratio[0])+0.5*(iTempN[0]-vData[i])*log(dZ2Ratio[0])
//				     -lgamma(0.5*(iTempN[0]+vData[i])+1)-lgamma(0.5*(iTempN[0]-vData[i])+1)-vLogSkellamZ[i]) <<endl;
//				cout << "CumU " << dCumU[0] << " dU "<< dU[0]<< endl;
//				cout<< "N "<<iTempN[0]<<endl;
//				cout<<"Skellam "<<vLogSkellamZ[i]<<endl;
//				cout<<"i "<<i<<"Data[i]"<<vData[i]<<endl;
//			}

		}

		sParam.vN[i]=iTempN[0];
		// cout<< "sParam.vN["<< i <<"] "<<sParam.vN[i]<<endl;
//		cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXx"<<endl;
		if(!isthisfinite( sParam.vN[i]))
		{
			cout <<"N enter"<<endl;
			double * vNError = new double [11];

			vNError[0]=i;
			vNError[1]=sParam.vN[i];
			vNError[2]=dCumU[0];
			vNError[3]=dU[0];
			vNError[4]=dLambda[0];
			vNError[5]=vLogSkellamZ[i];
			vNError[6]=vLogSkellamZ[i];
			vNError[7]=vData[i];
			vNError[8]=sParam.dMu[0];
			vNError[9]=sParam.vS[i];
			vNError[10]=sParam.vX[i];
			cout <<"N write"<<endl;
			string sFile="vNError.csv" ;
			WriteOutDoubleArray( vNError , 11, 1, sFile);
			cout <<"N del"<<endl;
			delete [] vNError;
			exit(1);

		}

	}

	delete dCumU;
	delete dLambda;
	delete iTempN;
	delete iNextBinom;
	delete dU;
	delete dZSum;
	delete dZ1Ratio;
	delete dZ2Ratio;
	delete dZLambda;


}

void DrawTauDNB(struct AllParam sParam,int iNumOfObs , const gsl_rng * gsl_random_num)
{

	double * dLambda=new double;
	for(int i=0; i<iNumOfObs; i++)
	{
//		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX "<< endl;
//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
//		cout << "sParam.vN[i] " << sParam.vN[i]<< endl;
//		cout<< "sParam.vX[i] "<< sParam.vX[i]<<endl;

		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);
		/* Drawing Beta(N,1)*/
		double dU2=gsl_rng_uniform(gsl_random_num);
		while(dU2==0)
		{
			dU2=gsl_rng_uniform(gsl_random_num);
		}
		if(sParam.vN[i]>0)
		{

			sParam.vTau2[i]=pow(dU2,(double) 1/sParam.vN[i]);
//			cout<< " 1/sParam.vN[i] "<< (double) 1/sParam.vN[i] << endl;
//			cout << "dU2 " << dU2<< endl;
//			cout << "sParam.vTau2[i] " << sParam.vTau2[i]<< endl;
		}
		else
		{
			sParam.vTau2[i]=0;
		}
		/* Drawing exponential with intensity 2 lambda*/

		double dU1=gsl_rng_uniform(gsl_random_num);
		sParam.vTau1[i]=1-log(1-dU1)/(dLambda[0]*(sParam.vZ1[i]+sParam.vZ2[i]))-sParam.vTau2[i];

		if(!isthisfinite(-log(sParam.vTau2[i])) & (sParam.vN[i]>0))
		{


			string sFile="vTau2Error.csv";
			WriteOutDoubleArray( sParam.vTau2, iNumOfObs, 1, sFile);
			sFile="dITau2Error.csv";
			WriteOutIntArray( &i, 1, 1, sFile);
			sFile="dUTau2Error.csv";
			WriteOutDoubleArray( &dU2, 1, 1, sFile);
			sFile="vNErrorTau.csv";
			WriteOutIntArray( sParam.vN, iNumOfObs, 1, sFile);
			sFile="vTau1Error.csv";
			WriteOutDoubleArray( sParam.vTau1, iNumOfObs, 1, sFile);
			exit(1);
		}


//		cout << "dU1 " << dU1<< endl;
//
//		cout << "sParam.vTau1[i] " << sParam.vTau1[i]<< endl;


	}
//	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX "<< endl;
}

void DrawAuxYandAuxHDNB(struct AllParam sParam, int iNumOfObs , const gsl_rng * gsl_random_num)
{
	int iNumOfComp;
	double* vWeights =new double [10];
	double* vMeans = new double [10];
	double* vVariances= new double [10];
	double* dU=new double;
	double* dCumP=new double[10];
	double* dLogLambda=new double;
	double* dLogZSum=new double;

	for(int k=0; k<iNumOfObs; k++)
	{

//		cout<<"XXXXXXXXXXXXXXXX  Start XXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;
//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
//		cout << "sParam.vN[i] " << sParam.vN[k]<< endl;
//		cout<< "sParam.vX[i] "<< sParam.vX[k]<<endl;
//		cout << "sParam.vTau1[k] " << sParam.vTau1[k]<<endl;
//		cout << "sParam.vTau2[k] " << sParam.vTau2[k]<<endl;

		dLogLambda[0]=sParam.dMu[0]+sParam.vS[k]+sParam.vX[k];
		dLogZSum[0]=log(sParam.vZ1[k]+sParam.vZ2[k]);
		/* Mixture approximation for tau1 */
		double mu = - digammaOnHost(1);
		double sigma = sqrt(trigammaOnHost(1));

		/* Weights */
		vWeights[0]=0.00396984425; vWeights[1]=0.0396244597; vWeights[2]=0.16776747; vWeights[3]=0.147036501; vWeights[4]=0.125306271;
		vWeights[5]=0.1014852; vWeights[6]=0.103758531; vWeights[7]=0.115972617; vWeights[8]=0.107065659; vWeights[9]=0.0880134482;

		/* Means */
		vMeans[0]=3.55887454*sigma+mu; vMeans[1]=2.11415904*sigma+mu; vMeans[2]= 0.968124631*sigma+mu; vMeans[3]= 0.51537638*sigma+mu; vMeans[4]= 0.145465449*sigma+mu;
		vMeans[5]=-0.145346445*sigma+mu; vMeans[6]= -0.416660312*sigma+mu; vMeans[7]=-0.689002855*sigma+mu; vMeans[8]= -0.974965634*sigma+mu; vMeans[9]= -1.27310004*sigma+mu;

		/* Variances */
		vVariances[0]=2.62603032*sigma*sigma; vVariances[1]=1.21263644*sigma*sigma; vVariances[2]=0.66586521*sigma*sigma; vVariances[3]=0.256650604*sigma*sigma; vVariances[4]= 0.120071142*sigma*sigma;
		vVariances[5]=0.0649909219*sigma*sigma; vVariances[6]=0.0473513798*sigma*sigma; vVariances[7]= 0.046639443*sigma*sigma; vVariances[8]= 0.0576541602*sigma*sigma; vVariances[9]= 0.0888536903*sigma*sigma;

//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//
//		}
//
//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//		}
//
//		for(int i=0 ; i<10; i++)
//		{
//			cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//		}

		/* Draw indicator form the mixture approximation for tau1*/
//		dCumP[0]=vWeights[0]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0],2) /vVariances[0])/sqrt(vVariances[0]);
		dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(-log(sParam.vTau1[k])-dLogZSum[0]-dLogLambda[0]-vMeans[0],2) /vVariances[0]);
//		dCumP[0]=exp(log(vWeights[0])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[0], sqrt(vVariances[0]))));
		for(int i=1; i<10;i++)
		{
//			dCumP[i]=dCumP[i-1]+vWeights[i]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i])/sqrt(vVariances[i]);
			dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(-log(sParam.vTau1[k])-dLogZSum[0]-dLogLambda[0]-vMeans[i],2) /vVariances[i]);
//			dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i], sqrt(vVariances[i]))));
//			cout<<"XX  contribution "<< vWeights[i]*exp(-0.5*pow(-log(sParam.vTau1[k])-log(2)-dLogLambda[0]-vMeans[i],2) /vVariances[i])/sqrt(vVariances[i])<<endl;
//			cout<<"XX  dCumP[i-1] "<< dCumP[i-1]<<endl;
//			cout<<"XX  dCumP[i] "<< dCumP[i]<<endl;

		}

//		for(int l=0; l<10; l++)
//		{
//			cout<<"1 dCumP[i]/dCumP[iNumOfComp-1] at "<<l << " is " << dCumP[l]/dCumP[9]<< endl;
//
//		}


		int i=0;
		dU[0]=gsl_rng_uniform(gsl_random_num);
//		cout<<"dU 1 " << dU[0]<<endl;
		while((dCumP[i]/dCumP[10-1])<dU[0])
		{
			i=i+1;
//			cout<<"i "<< i << endl;
		}

		/*Saving AuxY and AuxH for Tau1  */

//		cout<<"-log(sParam.vTau1[k]) "<<-log(sParam.vTau1[k])<<endl;
//		cout <<" -log 2 "<<-log(2)<<endl;
//		cout<<"-vMeans[i] "<<-vMeans[i] <<endl;
		sParam.mAuxY[2*k]=-log(sParam.vTau1[k])-dLogZSum[0]-vMeans[i];//-sParam.vS[k] ;
		sParam.mAuxH[4*k]=vVariances[i];


		if(sParam.vN[k]>0)
		{
			/* Mixture approximation for tau2 */
			MixtureLogGamma( sParam.vN[k], &iNumOfComp,  vWeights , vMeans ,vVariances);
			// cout<< "vN at "<<k<<" is " << sParam.vN[k] <<endl;
			// cout<< "iNumOfComp at "<<k<<" is " << iNumOfComp <<endl;

			/* Draw indicator form the mixture approximation for tau2*/
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//
//			}
//
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//			}
//
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//			}

//				dCumP[0]=vWeights[0]*exp(-0.5*pow((-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0]),2) /vVariances[0])/sqrt(vVariances[0]);
				dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(-log(sParam.vTau2[k])-dLogZSum[0]-dLogLambda[0]-vMeans[0],2) /vVariances[0]);
//				dCumP[0]=exp(log(vWeights[0])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[0], sqrt(vVariances[0]))));
				for(int i=1; i<iNumOfComp;i++)
				{
//					dCumP[i]=dCumP[i-1]+vWeights[i]*exp(-0.5*pow((-log(sParam.vTau2[k])-log( 2)-dLogLambda[0]-vMeans[i]),2) /vVariances[i])/sqrt(vVariances[i]);
					dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(-log(sParam.vTau2[k])-dLogZSum[0]-dLogLambda[0]-vMeans[i],2) /vVariances[i]);
//					dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])+log(gsl_ran_gaussian_pdf(-log(sParam.vTau2[k])-log(2)-dLogLambda[0]-vMeans[i], sqrt(vVariances[i]))));

//					cout<<"XX  contribution "<< vWeights[i]*exp(-0.5*pow((-log(sParam.vTau2[k])-log( 2)-dLogLambda[0]-vMeans[i]),2) /vVariances[i])/sqrt(vVariances[i])<<endl;
//					cout<<"XX  dCumP[i-1] "<< dCumP[i-1]<<endl;
//					cout<<"XX  dCumP[i] "<< dCumP[i]<<endl;
				}


//				for(int l=0; l<iNumOfComp; l++)
//				{
//					cout<<"sParam.vN[k] "<<sParam.vN[k]<<" iNumOfComp "<< iNumOfComp << endl;
//					cout<<"2 dCumP[i]/dCumP[iNumOfComp-1] at "<<l << " is " << dCumP[l]/dCumP[iNumOfComp-1]<< endl;
//
//				}

				int j=0;
				dU[0]=gsl_rng_uniform(gsl_random_num);
//				cout<<"dU 2 " << dU[0]<<endl;
				while(dCumP[j]/dCumP[iNumOfComp-1]<dU[0])
				{
					j=j+1;
//					cout<<"j "<< j << endl;
				}

				/*Saving AuxY and AuxH for Tau2 */

				if(!isthisfinite( -log(sParam.vTau2[k])-dLogZSum[0]-vMeans[j]))
										{
											double * vYErrorAux = new double [4];

											vYErrorAux[0]=k;
											vYErrorAux[1]=sParam.vTau2[k];
											vYErrorAux[2]=-log(sParam.vTau2[k]);
											vYErrorAux[3]=-dLogZSum[0];

											cout <<" sParam.vTau2[k] "<< sParam.vTau2[k]<<endl;
											cout <<" -log(sParam.vTau2[k]) "<< -log(sParam.vTau2[k])<<endl;
											cout <<" -dLogZSum[0] "<< -dLogZSum[0]<<endl;
											cout<<" z1 "<<sParam.vZ1[k] <<endl;
											cout<<" z2 "<<sParam.vZ2[k] <<endl;
											cout <<" -vMeans[j] "<< -vMeans[j]<<endl;

											string sFile="vYErrorAux.csv";
											WriteOutDoubleArray( vYErrorAux , 4, 1, sFile);

											delete [] vYErrorAux;
											cout<<"AuxY Error"<<endl;
											exit(1);

										}


				if(!isthisfinite( 1.0/vVariances[j]))
								{
									double * vVarErrorAux = new double [iNumOfComp+2];

									for(int t=0; t<iNumOfComp;t++)
									{
										vVarErrorAux[t]=vVariances[t];
									}

									vVarErrorAux[iNumOfComp]=sParam.vN[k];
									vVarErrorAux[iNumOfComp+1]=k;
									string sFile="vVarErrorAux.csv";
									WriteOutDoubleArray( vVarErrorAux , iNumOfComp+2, 1, sFile);

									cout<<"AuxH Error"<<endl;
									delete [] vVarErrorAux;
									exit(1);
					}


				sParam.mAuxY[2*k+1]= -log(sParam.vTau2[k])-dLogZSum[0]-vMeans[j];//-sParam.vS[k];
				sParam.mAuxH[4*k+3]=vVariances[j] ;







		}
		else
		{
			sParam.mAuxY[2*k+1]= 0;
			sParam.mAuxH[4*k+3]=0;

		}


//		cout<<"AuxY 0 at "<<k << " is " <<sParam.mAuxY[2*k] << endl;
//		cout<<"AuxY 1 at "<<k << " is " <<sParam.mAuxY[2*k+1] << endl;
//		cout<<"AuxH 0 at "<<k << " is " <<sParam.mAuxH[4*k] << endl;
//		cout<<"AuxH 1 at "<<k << " is " <<sParam.mAuxH[4*k+3] << endl;
//
//		cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;

		sParam.mAuxH[4*k+1]=0;
		sParam.mAuxH[4*k+2]=0;


	}

	delete [] vWeights;
	delete [] vMeans;
	delete [] vVariances;
	delete [] dCumP;
	delete dLogZSum;
	delete dU;
	delete dLogLambda;

}

void GammaLogPosteriorDNB(double * dLogGamma, double *dPriorA, double *dPriorB,  int iNumOfObs,
		double * vDNB, double* vIndicator, double * dLogPosterior)
{
	dLogPosterior[0]=0;
	double dGamma=LogitTransformBack(dLogGamma[0]);

	for(int i=0; i<iNumOfObs; i++)
	{
		dLogPosterior[0]+=log( dGamma*vIndicator[i]+(1-dGamma)*vDNB[i] );
	}
	dLogPosterior[0]+=(dPriorA[0] -1)*log(dGamma)+(dPriorB[0] -1)*log(1-dGamma);
}

void DrawGammaAdaptiveRWDNB(  int iNumOfObs,  int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam,
		double* vIndicator, double* vDNB, double* dCovar,  double* dSum  )
{


	/* Proposal */
	double dOmega1=0.05;
	double dU=gsl_rng_uniform(gsl_random_num);
	double dNewGammaLogit;
	double dSigma;
	double * dLogPosteriorNew=new double;
	double * dLogPosteriorOld=new double;

	if(iNumOfIter<=1)
	{
		dSigma=0.1;
	}
	else if(iNumOfIter==2)
	{
		dSigma=pow(dCovar[0]-0.5*dSum[0]*dSum[0],0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}
	else
	{
		dSigma=pow(dCovar[0], 0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}

	if(dU<=dOmega1)
	{
		/* Static proposal */
		dNewGammaLogit=LogitTransform(sParam.dGamma[0])+gsl_ran_gaussian(gsl_random_num,0.1);
	}
	else
	{
		/* Proposal based on empirical covariance */
		dNewGammaLogit=LogitTransform(sParam.dGamma[0])+gsl_ran_gaussian(gsl_random_num,dSigma*2.38);
	}


	/* Acceptance Rate */

	double dOldGammaLogit=LogitTransform(sParam.dGamma[0]);

	GammaLogPosteriorDNB(&dNewGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vDNB, vIndicator, dLogPosteriorNew)	;
	GammaLogPosteriorDNB(&dOldGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vDNB, vIndicator, dLogPosteriorOld)	;

	double dLogU=log(gsl_rng_uniform(gsl_random_num));
    double dAccept=fmin(dLogPosteriorNew[0]-dLogPosteriorOld[0],0);

	if((dLogU<=dAccept) & (isthisfinite(LogitTransformBack(dNewGammaLogit))) )
	{
		sParam.dGamma[0]=LogitTransformBack(dNewGammaLogit);
	}


	/* Covariance */
	if(iNumOfIter<=1)
	{
		dCovar[0]=dCovar[0]+pow(LogitTransform(sParam.dGamma[0]),2);
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
	}
	else if(iNumOfIter==2)
	{
		dCovar[0]=(iNumOfIter-1)*pow(dSigma,2)/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));
	}
	else
	{
		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+(dSum[0]/iNumOfIter)*(dSum[0]/iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-(dSum[0]/iNumOfIter)*(dSum[0]/(iNumOfIter+1));
	}

	cout << ":: ARW2 dSum " << dSum[0]<< endl;
	cout << ":: ARW2 dSigma " << dSigma<< endl;
	cout << ":: ARW2 dCovar " << dCovar[0]<< endl;
	cout << ":: ARW proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
	cout << ":: ARW accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;

}

void NuLogPosteriorDNB(double * dLogNu, double *vZ1, double* vZ2, double *dPriorA, double *dPriorB,  int iNumOfObs,
	double * dLogPosterior)
{

	double dNu=exp(dLogNu[0]);
	dLogPosterior[0]=2*iNumOfObs*(dNu*log(dNu)-lgamma(dNu) );

	for(int i=0; i<iNumOfObs; i++)
	{
		dLogPosterior[0]+=dNu*(log(vZ1[i])+ log(vZ2[i]) -vZ1[i]-vZ2[i] );
	}
	dLogPosterior[0]+=dPriorA[0]*log(dPriorB[0]) +dPriorA[0]*log(dNu) -dPriorB[0]*dNu-lgamma(dPriorA[0]);


}

void DrawNuAdaptiveRWDNB( int iNumOfObs,  int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam, double* dCovar,  double* dSum  )
{


	/* Proposal */
	double dOmega1=0.05;
	double dU=gsl_rng_uniform(gsl_random_num);
	double dNewNuLog;
	double dSigma;
	double * dLogPosteriorNew=new double;
	double * dLogPosteriorOld=new double;

	if(iNumOfIter<=1)
	{
		dSigma=0.1;
	}
	else if(iNumOfIter==2)
	{
		dSigma=pow(dCovar[0]-0.5*dSum[0]*dSum[0],0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}
	else
	{
		dSigma=pow(dCovar[0], 0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}

	if(dU<=dOmega1)
	{
		/* Static proposal */
//		dNewNuLog=log(sParam.dNu[0])+gsl_ran_gaussian(gsl_random_num,0.1);
		dNewNuLog=log(sParam.dNu[0])+gsl_ran_tdist(gsl_random_num,5)*0.1;
	}
	else
	{
		/* Proposal based on empirical covariance */
//		dNewNuLog=log(sParam.dNu[0])+gsl_ran_gaussian(gsl_random_num,dSigma*2.38);
		dNewNuLog=log(sParam.dNu[0])+gsl_ran_tdist(gsl_random_num,5)*dSigma*2.38;
	}


	/* Acceptance Rate */

	double dOldNuLog=log(sParam.dNu[0]);

	NuLogPosteriorDNB(&dNewNuLog,sParam.vZ1,sParam.vZ2, sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorNew)	;
	NuLogPosteriorDNB(&dOldNuLog, sParam.vZ1,sParam.vZ2, sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorOld)	;

//	cout <<"dLogPosteriorNew[0] "<<dLogPosteriorNew[0]<<endl;
//	cout <<"dLogPosteriorOld[0] "<<dLogPosteriorOld[0]<<endl;

	double dLogU=log(gsl_rng_uniform(gsl_random_num));
    double dAccept=fmin(dLogPosteriorNew[0]-dLogPosteriorOld[0],0);

//    cout <<"dAccpet "<<dAccept<<" dLogU "<< dLogU<<endl;
	if(dLogU<=dAccept)
	{
		sParam.dNu[0]=exp(dNewNuLog);
//		cout <<"Accepted! "<<endl;
	}



	/* Covariance */
	if(iNumOfIter<=1)
	{
		dCovar[0]=dCovar[0]+pow(log(sParam.dNu[0]),2);
		dSum[0]=dSum[0]+log(sParam.dNu[0]);

	}
	else if(iNumOfIter==2)
	{
		dCovar[0]=(iNumOfIter-1)*pow(dSigma,2)/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(log(sParam.dNu[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+log(sParam.dNu[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));


	}
	else
	{
		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(log(sParam.dNu[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+log(sParam.dNu[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));
	}

//	cout << "dSum " << dSum[0]<< endl;
//	cout << "dCovar " << dCovar[0]<< endl;
//	cout << "proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
//	cout << "accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;

}

double NuLogPosteriorGSL(const gsl_vector *v, void *ParamPointer)
{
	double dValue=0;
	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double dLogNu= gsl_vector_get(v, 0);

	double dNu=exp(dLogNu);
	dValue=2*(Param->iNumOfObs)[0]*(dNu*log(dNu)-lgamma(dNu) );

	for(int i=0; i<(Param->iNumOfObs)[0]; i++)
	{
		dValue+=dNu*(log((Param->vZ1)[i])+ log((Param->vZ2)[i]) -(Param->vZ1)[i]-(Param->vZ2)[i] );
	}
	dValue+=(Param->dPriorNuA)[0]*log((Param->dPriorNuB)[0]) +(Param->dPriorNuA)[0]*log(dNu) -(Param->dPriorNuB)[0]*dNu-lgamma((Param->dPriorNuA)[0]);

	return -dValue/((Param->iNumOfObs)[0]);
}

void NuLogPosteriorDerivativeGSL(const gsl_vector *v, void *ParamPointer,
       gsl_vector *df)
{
	double dValue=0;
	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double dLogNu= gsl_vector_get(v, 0);
	double dNu=exp(dLogNu);
	dValue=2*(Param->iNumOfObs)[0]*(log(dNu)+1-gsl_sf_psi(dNu) );

	for(int i=0; i<(Param->iNumOfObs)[0]; i++)
	{
		dValue+=(log((Param->vZ1)[i])+ log((Param->vZ2)[i]) -(Param->vZ1)[i]-(Param->vZ2)[i] );
	}
	dValue+=(Param->dPriorNuA)[0]/dNu -(Param->dPriorNuB)[0];

	gsl_vector_set(df, 0, -dValue/((Param->iNumOfObs)[0]));
}

void NuLogPosteriorANDDerivativeGSL(const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
  *f = NuLogPosteriorGSL(x, params);
  NuLogPosteriorDerivativeGSL(x, params, df);
}

void DrawNuLaplaceApproxDNB( struct AllParam sParam , int iNumOfObs, const gsl_rng * gsl_random_num)
{
	/* Setting up the optimizer */
	int iter = 0;
	int max_iter=100;
	int status;

//	cout << "Num Of obs : " << iNumOfObs << endl;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf F;


//	cout << "Num Of obs 2: " << iNumOfObs << endl;
	/* setting up the objective function */
	F.n = 1;
	F.f =&NuLogPosteriorGSL;
	F.df = &NuLogPosteriorDerivativeGSL;
	F.fdf = &NuLogPosteriorANDDerivativeGSL;
	F.params = &sParam;

	/* Starting point*/
	x = gsl_vector_alloc (1);
	gsl_vector_set (x,0, log(5));

	/* Setting up the minimizer */
	T =  gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, 1);

	gsl_multimin_fdfminimizer_set (s, &F, x, 0.000001, 0.000001);

	/* Iterating the optimizer */
	do
	{
	      iter++;
	      status = gsl_multimin_fdfminimizer_iterate (s);

//	      cout<<iter<<endl;
	      if (status)
	        break;

	      status = gsl_multimin_test_gradient (s->gradient, 1e-2);

//	      if (status == GSL_SUCCESS)
//	      printf ("Minimum found at:\n");
//
//	      printf ("%5d %.5f %10.5f\n", iter,
//	              gsl_vector_get (s->x, 0),
//	              s->f);

	  }
	  while (status == GSL_CONTINUE && iter < max_iter);


//	  /* Mean */
	  double dMuLogit = gsl_vector_get (s->x, 0);


	  gsl_vector *dfx;
	  gsl_vector *dfxh;
	  gsl_vector *xh;
	  dfx=gsl_vector_alloc (1);
	  dfxh=gsl_vector_alloc (1);
	  xh=gsl_vector_alloc (1);
	  double dTol=0.0000001;

	  //cout <<  "Logit dMu:  "<<gsl_vector_get (s->x, 0) << endl;
	  //cout << "dMu: "<<LogitTransformBack(gsl_vector_get (s->x, 0)) << endl;

	  //cout <<  "Logit dMu +h: "<<gsl_vector_get (s->x, 0)+dTol << endl;
	  //cout << "dMu +h: "<<LogitTransformBack(gsl_vector_get (s->x, 0)+dTol) << endl;

	  gsl_vector_set (xh,0, gsl_vector_get (s->x, 0)+dTol);


	  NuLogPosteriorDerivativeGSL(s->x,  &sParam ,dfx);
	  NuLogPosteriorDerivativeGSL(xh,  &sParam ,dfxh);

	  /* Sigma */

	  double dSigmaLogit=sqrt(dTol/((gsl_vector_get (dfxh, 0)-gsl_vector_get (dfx, 0))*(iNumOfObs)));

	  //cout<< "dfxh: " << gsl_vector_get (dfxh,0) << endl;
	  //cout<< "dfx: " << gsl_vector_get (dfx,0)<< endl;
	  //cout<< "dSigmaLog: " << dSigmaLogit << endl;
	  //cout<< "dMu: " << dMuLogit << endl;
	  /* Draw new Gamma log */
	  //double dNewGammaLogit=dMuLogit+gsl_ran_gaussian(gsl_random_num,dSigmaLogit);
	  double dNewNuLogit=dMuLogit+gsl_ran_tdist(gsl_random_num,4)*dSigmaLogit;

	  gsl_vector * newx;
	  newx=gsl_vector_alloc(1);
	  gsl_vector_set (newx,0, dNewNuLogit);

	  gsl_vector * oldx;
	  oldx=gsl_vector_alloc(1);
	  gsl_vector_set (oldx,0, log(sParam.dNu[0]));




	  /*Calculate Acceptance rate */
	  double dLogU=log(gsl_rng_uniform(gsl_random_num));
	  double dAccept=fmin(-NuLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(gsl_ran_tdist_pdf((log(sParam.dNu[0])-dMuLogit)/dSigmaLogit,4))
					 		 +NuLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(gsl_ran_tdist_pdf((dNewNuLogit-dMuLogit)/dSigmaLogit,4)),0);

	  if(dLogU<=dAccept)
	  {
		  sParam.dNu[0]=exp(dNewNuLogit);

	  }

	  gsl_multimin_fdfminimizer_free (s);
	  gsl_vector_free (x);



}

void DrawNuDiscreteDNB(struct AllParam sParam ,  int iNumOfObs, const gsl_rng * gsl_random_num, double dResolution,  double dDelta)
{
	/* Propose draw */
	double dU=gsl_rng_uniform(gsl_random_num);
	int iN=2*dDelta/dResolution+1;

	double dCum=1.0/iN;
	double dNuProposed=sParam.dNu[0]-dDelta;
	while(dCum<dU)
	{
		dCum=dCum+1.0/iN;
		dNuProposed=dNuProposed+dResolution;
	}

//	cout<< "dNuProposed "<<dNuProposed<<endl;
//	cout<< "dNuCurrent "<<sParam.dNu[0]<<endl;

	if((dNuProposed < 128) & (dNuProposed >2))
	{
		/* Calculate acceptance probability */
		double dLogProposal=2*iNumOfObs*(dNuProposed*log(dNuProposed)-lgamma(dNuProposed) );
		double dLogCurrent=2*iNumOfObs*(sParam.dNu[0]*log(sParam.dNu[0])-lgamma(sParam.dNu[0]) );
		for(int i=0; i<iNumOfObs; i++)
		{
			dLogProposal=dLogProposal+dNuProposed*(log(sParam.vZ1[i])+ log(sParam.vZ2[i]) -sParam.vZ1[i]-sParam.vZ2[i]);
			dLogCurrent=dLogCurrent+sParam.dNu[0]*(log(sParam.vZ1[i])+ log(sParam.vZ2[i]) -sParam.vZ1[i]-sParam.vZ2[i]);
		}

//		cout<< "dLogProposed "<<dLogProposal<<endl;
//		cout<< "dLogCurrent "<<dLogCurrent<<endl;

		double dAccept=fmin(dLogProposal-dLogCurrent,0);

		/* Accept or reject*/
		double dLogU=log(gsl_rng_uniform(gsl_random_num));

//		cout << "dAccept "<<dAccept<<" dLogU "<< dLogU<<endl;
		if(dLogU<=dAccept)
		{
			sParam.dNu[0]=dNuProposed;
		}
	}

}

void CalculateWStar( int iNumOfObs, double* vTimes,  int iNumOfKnots, double* vKnots,gsl_matrix * mLambdaInvTheta,  gsl_matrix * mWStar )
{

	double *dH=new double;

	list<double> vUniqueTimes;
	for(int i=0; i<iNumOfObs; i++)
	{
		vUniqueTimes.push_back(vTimes[i]);
	}


	vUniqueTimes.sort();
	vUniqueTimes.unique();



//	for (std::list<double>::iterator it=vUniqueTimes.begin(); it!=vUniqueTimes.end(); ++it)
//	{
//		  std::cout <<  *it << endl;
//	}

	int iNumOfUnique=vUniqueTimes.size();

	gsl_matrix * mPUnique = gsl_matrix_alloc (iNumOfUnique, iNumOfKnots);
	gsl_matrix * mWUnique  = gsl_matrix_alloc (iNumOfUnique, iNumOfKnots);

	list<double>::iterator it=vUniqueTimes.begin();
	for(int i=0; i<iNumOfUnique;i++)
	{

		int iTempKnot=1;
		while(vTimes[i] > vKnots[iTempKnot])
		{

			iTempKnot=iTempKnot+1;
		}

		for(int j=0; j<iNumOfKnots; j++)
		{
			if(j==iTempKnot-1)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mPUnique , i, j, (vKnots[iTempKnot]- *it )*(pow(vKnots[iTempKnot]-*it ,2)-dH[0]*dH[0])/(6*dH[0]) );
				gsl_matrix_set (mWUnique , i, j, (vKnots[iTempKnot]- *it )/dH[0] );

			}
			else if (j==iTempKnot)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mPUnique , i, j, (*it -vKnots[iTempKnot-1])*(pow( *it -vKnots[iTempKnot-1],2)-dH[0]*dH[0])/(6*dH[0]));
				gsl_matrix_set (mWUnique , i, j, (*it -vKnots[iTempKnot-1])/dH[0]);

			}
			else
			{
				gsl_matrix_set (mPUnique , i, j, 0);
				gsl_matrix_set (mWUnique , i, j, 0);
			}
		}

		it++;
	}

//	for(int i=0; i<iNumOfObs; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			 printf ("mPUnique(%d,%d) = %g\n", i, j, gsl_matrix_get (mPUnique, i, j));
//		}
//	}

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mPUnique,mLambdaInvTheta,1.0, mWUnique);


	for(int j=0; j<iNumOfKnots; j++)
	{
		gsl_matrix_set (mWStar ,0,j,0);
		for(int i=0; i<iNumOfUnique; i++)
		{
			gsl_matrix_set (mWStar ,0,j, gsl_matrix_get (mWStar ,0,j)+gsl_matrix_get (mWUnique ,i,j) );
		}
	}


	delete dH;
}

void CalculateSplineMatrices( int iNumOfObs, double* vTimes,  int iNumOfKnots, double* vKnots, gsl_matrix * mWTilde)
{
	gsl_matrix * mTheta = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_matrix * mLambda = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_matrix * mP = gsl_matrix_alloc (iNumOfObs, iNumOfKnots);
	gsl_matrix * mW = gsl_matrix_alloc (iNumOfObs, iNumOfKnots);
	gsl_matrix * mWStar = gsl_matrix_alloc (1, iNumOfKnots);



	double *dLambda=new double;
	double *dH=new double;
	double *dHPlus=new double;
	/* Lambda and Theta */
	for(int i=0; i<iNumOfKnots; i++)
	{
		for(int j=0; j<iNumOfKnots; j++)
		{
			gsl_matrix_set (mLambda, i, j, 0);
			gsl_matrix_set (mTheta, i, j, 0);
		}

		if(i==0 ||  i==iNumOfKnots-1 )
		{
			gsl_matrix_set (mLambda , i, i, 2);

		}
		else
		{
			dH[0]=vKnots[i]-vKnots[i-1];
			dHPlus[0]=vKnots[i+1]-vKnots[i];
			dLambda[0]=dHPlus[0]/(dH[0]+dHPlus[0]);
			/* Lambda matrix */
			gsl_matrix_set (mLambda , i, i-1, 1-dLambda[0]);
			gsl_matrix_set (mLambda, i, i, 2);
			gsl_matrix_set (mLambda , i, i+1, dLambda[0]);
			/* Theta matrix */
			gsl_matrix_set (mTheta , i, i-1, 6/(dH[0]*(dH[0]+dHPlus[0])));
			gsl_matrix_set (mTheta , i, i, -6/(dH[0]*dHPlus[0]));
			gsl_matrix_set (mTheta , i, i+1,6/(dHPlus[0]*(dH[0]+dHPlus[0])));

		}

	}

	/* P and Q */
	for(int i=0; i<iNumOfObs;i++)
	{
		int iTempKnot=1;
		while(vTimes[i] > vKnots[iTempKnot])
		{
//			cout<<"vTimes "<<vTimes[i]<<"vKnots "<<vKnots[iTempKnot]<<endl;
			iTempKnot=iTempKnot+1;
		}
//		cout <<"END "<<endl;
		for(int j=0; j<iNumOfKnots; j++)
		{
			if(j==iTempKnot-1)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mP , i, j, (vKnots[iTempKnot]-vTimes[i])*(pow(vKnots[iTempKnot]-vTimes[i],2)-dH[0]*dH[0])/(6*dH[0]) );
				gsl_matrix_set (mW , i, j, (vKnots[iTempKnot]-vTimes[i])/dH[0] );

			}
			else if (j==iTempKnot)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mP , i, j, (vTimes[i]-vKnots[iTempKnot-1])*(pow(vTimes[i]-vKnots[iTempKnot-1],2)-dH[0]*dH[0])/(6*dH[0]));
				gsl_matrix_set (mW , i, j, (vTimes[i]-vKnots[iTempKnot-1])/dH[0]);

			}
			else
			{
				gsl_matrix_set (mP , i, j, 0);
				gsl_matrix_set (mW , i, j, 0);
			}
		}
	}

//	for(int i=0; i<iNumOfKnots; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			 printf ("mLambda(%d,%d) = %g\n", i, j, gsl_matrix_get (mLambda, i, j));
//
//		}
//	}
//
//	for(int i=0; i<iNumOfKnots; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			 printf ("mTheta(%d,%d) = %g\n", i, j, gsl_matrix_get (mTheta, i, j));
//
//		}
//	}
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			 printf ("mP(%d,%d) = %g\n", i, j, gsl_matrix_get (mP, i, j));
//		}
//	}
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			 printf ("mQ(%d,%d) = %g\n", i, j, gsl_matrix_get (mW, i, j));
//		}
//	}
//
//	for(int i=4764; i<4768;  i++)
//	{
//
//		printf ("vTimes(%d) = %g\n", i, vTimes[i]);
//
//	}
//
//	for(int i=4764; i<4768; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			printf ("mQ(%d,%d) = %g\n", i, j, gsl_matrix_get (mW, i, j));
//		}
//	}
//
//	for(int i=4764; i<4768;  i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			printf ("mP(%d,%d) = %g\n", i, j, gsl_matrix_get (mP, i, j));
//		}
//	}

	/* Lambda  inverse */
	gsl_matrix * mLambdaInv = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots);

	int s;
	gsl_linalg_LU_decomp (mLambda, mPermutation, &s);
	gsl_linalg_LU_invert (mLambda, mPermutation, mLambdaInv);



	/* W matrix */
	gsl_matrix * mLambdaInvTheta = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	for(int i=0; i<iNumOfKnots; i++)
	{
		for(int j=0; j<iNumOfKnots; j++)
		{
			gsl_matrix_set ( mLambdaInvTheta  ,1,j,0 );
		}
	}

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mLambdaInv,mTheta,0.0, mLambdaInvTheta);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mP,mLambdaInvTheta,1.0, mW);

//	/* Old WStar*/
//	gsl_matrix * mWStarOld = gsl_matrix_alloc (1, iNumOfKnots);
//	for(int j=0; j<iNumOfKnots; j++)
//	{
//		gsl_matrix_set (mWStarOld ,0,j,0);
//		for(int i=0; i<iNumOfObs; i++)
//		{
//			gsl_matrix_set (mWStarOld ,0,j, gsl_matrix_get (mWStarOld,0,j)+gsl_matrix_get (mW ,i,j) );
//		}
//	}


	/* WStar matrix */
//	cout <<"calculate star started "<<endl;
	CalculateWStar( iNumOfObs, vTimes, iNumOfKnots,  vKnots, mLambdaInvTheta, mWStar );
//	cout <<"calculate star finished "<<endl;


	/* WTilde matrix */
	for(int i=0; i<iNumOfObs; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			gsl_matrix_set (mWTilde ,i,j, gsl_matrix_get (mW ,i,j) - gsl_matrix_get (mW ,i,iNumOfKnots-1)*gsl_matrix_get (mWStar ,0,j)/gsl_matrix_get (mWStar ,0,iNumOfKnots-1) );
		}
	}


	/* PRINT PRINT */
//	for(int i=0; i<iNumOfKnots; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//
//			 printf ("mLambda Inv(%d,%d) = %g\n", i, j, gsl_matrix_get (mLambdaInv, i, j));
//		}
//	}
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			 printf ("mW(%d,%d) = %g\n", i, j,
//				              gsl_matrix_get (mW, i, j));
//		}
//	}
//
//	for(int j=0; j<iNumOfKnots; j++)
//	{
//
//		cout<<"WStat "<<std::setprecision(12)<<  gsl_matrix_get (mWStar, 0, j)<<endl;
//	}
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		for(int j=0; j<iNumOfKnots-1; j++)
//		{
//			 printf ("mWTilde(%d,%d) = %g\n", i, j,
//			              gsl_matrix_get (mWTilde, i, j));
//		}
//	}
//
//	for(int i=4764; i<4768;  i++)
//	{
//		for(int j=0; j<iNumOfKnots; j++)
//		{
//			printf ("mW(%d,%d) = %g\n", i, j, gsl_matrix_get (mW, i, j));
//		}
//	}
//
//	for(int i=4764; i<4768;  i++)
//	{
//		for(int j=0; j<iNumOfKnots-1; j++)
//		{
//			printf ("mWtilde(%d,%d) = %g\n", i, j, gsl_matrix_get (mWTilde, i, j));
//		}
//	}

	gsl_matrix_free (mTheta);
	gsl_matrix_free (mLambda);
	gsl_matrix_free (mP);
	gsl_matrix_free (mW);
	gsl_matrix_free (mWStar);
	gsl_matrix_free (mLambdaInv);

	gsl_permutation_free (mPermutation);

	delete dLambda;
	delete dH;
	delete dHPlus;

}

void DrawSesonal(struct AllParam sParam,  int iNumOfObs, int iNumOfKnots,   gsl_matrix * mWTilde ,const gsl_rng * gsl_random_num)
{

	/* Number of observations in the OLS */
	int iNumOfOLSObs=0;
	for(int i=0; i<iNumOfObs;i++)
	{
		if(sParam.vN[i]==0)
		{
			iNumOfOLSObs=iNumOfOLSObs+1;
		}
		else
		{
			iNumOfOLSObs=iNumOfOLSObs+2;
		}
	}
	gsl_matrix * mYBar=gsl_matrix_alloc (iNumOfOLSObs, 1);
	gsl_matrix * mWBar=gsl_matrix_alloc(iNumOfOLSObs,iNumOfKnots-1);

	double * dY1=new double;
	double * dY2=new double;

	double * dS1=new double;
	double * dS2=new double;


	int iTempCount=0;
	for(int i=0; i<iNumOfObs; i++)
	{


		if(sParam.vN[i]==0)
		{

			dY1[0]=sParam.mAuxY[i*2]-(sParam.dMu[0]+sParam.vX[i]);
			dS1[0]=pow(sParam.mAuxH[i*4],0.5);

			gsl_matrix_set (mYBar,iTempCount,0, dY1[0]/dS1[0] );

			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,iTempCount,j, gsl_matrix_get (mWTilde, i, j)/dS1[0]);
			}

			iTempCount=iTempCount+1;


		}
		else
		{

			dY1[0]=sParam.mAuxY[i*2]-(sParam.dMu[0]+sParam.vX[i]);
			dY2[0]=sParam.mAuxY[i*2+1]-(sParam.dMu[0]+sParam.vX[i]);
			dS1[0]=pow(sParam.mAuxH[i*4],0.5);
			dS2[0]=pow(sParam.mAuxH[i*4+3],0.5);

			gsl_matrix_set (mYBar,iTempCount,0, dY1[0]/dS1[0]);
			gsl_matrix_set (mYBar,iTempCount+1,0, dY2[0]/dS2[0]);

			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,iTempCount,j, gsl_matrix_get (mWTilde, i, j)/dS1[0]);
				gsl_matrix_set (mWBar,iTempCount+1,j, gsl_matrix_get (mWTilde, i, j)/dS2[0]);
			}

			iTempCount=iTempCount+2;
		}

	}




	/* Variance */
	gsl_matrix * mWW=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVar=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVarInv=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mWBar,0.0, mWW);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			if(j==i)
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j)+sParam.dPriorSeasonalVar[0]);
			}
			else
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j));
			}
		}
	}

	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots-1);
	int s;
	gsl_linalg_LU_decomp (mVarInv, mPermutation, &s);
	gsl_linalg_LU_invert (mVarInv, mPermutation, mVar);



	/* Mean */
	gsl_matrix * mWY=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMean=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMeanTemp=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mYBar,0.0, mWY);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		gsl_matrix_set (mMeanTemp,i,0,gsl_matrix_get(mWY,i,0)+sParam.dPriorSeasonalMean[0]/sParam.dPriorSeasonalVar[0]);
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mVar,mMeanTemp,0.0, mMean);


	/* Draw Normal */
	gsl_linalg_cholesky_decomp (mVar);
	gsl_matrix * vNormal=gsl_matrix_alloc (iNumOfKnots-1, 1);


	string sFile;
	for(int i=0; i<iNumOfKnots-1;i++)
	{

		if(!isthisfinite(gsl_matrix_get(mMean,i,0)))
		{
			double * YBar=new double[iNumOfOLSObs];
			for(int k=0; k<iNumOfOLSObs;k++)
			{
				YBar[k]=gsl_matrix_get(mYBar,k,0);
			}
			sFile="vYBarSpineError.csv";
			WriteOutDoubleArray( YBar , iNumOfOLSObs, 1, sFile);

			double * WBar=new double[iNumOfOLSObs*(iNumOfKnots-1)];
			for(int k=0; k<iNumOfOLSObs;k++)
			{
				for(int l=0; l<iNumOfKnots-1;l++)
				{
					WBar[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mWBar,k,l);
				}
			}
			sFile="vWBarSpineError.csv";
			WriteOutDoubleArray( WBar , iNumOfOLSObs, iNumOfKnots-1, sFile);

			double * Mean=new double[iNumOfKnots-1];
			for(int k=0; k<iNumOfKnots-1;k++)
			{
				Mean[k]=gsl_matrix_get(mMean,k,0);
			}
			sFile="vMeanSpineError.csv";
			WriteOutDoubleArray( Mean , iNumOfKnots-1, 1, sFile);





			sFile="vHSpineError.csv";
			WriteOutDoubleArray( sParam.mAuxH , iNumOfObs, 4, sFile);


			sFile="vAuxYSpineError.csv";
			WriteOutDoubleArray( sParam.mAuxY , iNumOfObs, 2, sFile);

			sFile="vXSpineError.csv";
			WriteOutDoubleArray( sParam.vX , iNumOfObs, 1, sFile);

			sFile="vSSpineError.csv";
			WriteOutDoubleArray( sParam.vS , iNumOfObs, 1, sFile);

			double * Var=new double[(iNumOfKnots-1)*(iNumOfKnots-1)];
			for(int k=0; k<iNumOfKnots-1;k++)
			{
				for(int l=0; l<iNumOfKnots-1;l++)
				{
					Var[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mVar,k,l);
				}
			}
			sFile="vVarSpineError.csv";
			WriteOutDoubleArray( Var , iNumOfKnots-1, iNumOfKnots-1, sFile);

			delete [] YBar;
			delete [] WBar;
			delete [] Mean;
			delete [] Var;

			exit(1);
		}

		for(int j=0; j<iNumOfKnots-1;j++)
		{
			if(!isthisfinite(gsl_matrix_get(mVar,i,j)))
			{
				double * YBar=new double[iNumOfOLSObs];
				for(int k=0; k<iNumOfOLSObs;k++)
				{
					YBar[k]=gsl_matrix_get(mYBar,k,0);
				}
				sFile="vYBarSpineError.csv";
				WriteOutDoubleArray( YBar , iNumOfOLSObs, 1, sFile);

				double * WBar=new double[iNumOfOLSObs*(iNumOfKnots-1)];
				for(int k=0; k<iNumOfOLSObs;k++)
				{
					for(int l=0; l<iNumOfKnots-1;l++)
					{
						WBar[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mWBar,k,l);
					}
				}
				sFile="vWBarSpineError.csv";
				WriteOutDoubleArray( WBar , iNumOfOLSObs, iNumOfKnots-1, sFile);

				double * Mean=new double[iNumOfKnots-1];
				for(int k=0; k<iNumOfKnots-1;k++)
				{
					Mean[k]=gsl_matrix_get(mMean,k,0);
				}
				sFile="vMeanSpineError.csv";
				WriteOutDoubleArray( Mean , iNumOfKnots-1, 1, sFile);

				sFile="vHSpineError.csv";
				WriteOutDoubleArray( sParam.mAuxH , iNumOfObs, 4, sFile);

				sFile="vAuxYSpineError.csv";
				WriteOutDoubleArray( sParam.mAuxY , iNumOfObs, 2, sFile);

				sFile="vXSpineError.csv";
				WriteOutDoubleArray( sParam.vX , iNumOfObs, 1, sFile);

				sFile="vSSpineError.csv";
				WriteOutDoubleArray( sParam.vS , iNumOfObs, 1, sFile);

				double * Var=new double[(iNumOfKnots-1)*(iNumOfKnots-1)];
				for(int k=0; k<iNumOfKnots-1;k++)
				{
					for(int l=0; l<iNumOfKnots-1;l++)
					{
						Var[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mVar,k,l);
					}
				}
				sFile="vVarSpineError.csv";
				WriteOutDoubleArray( Var , iNumOfKnots-1, iNumOfKnots-1, sFile);

				delete [] YBar;
				delete [] WBar;
				delete [] Mean;
				delete [] Var;
				exit(1);
			}
		}

	}

	/* Draw Beta */
	gsl_matrix * vBeta=gsl_matrix_alloc (iNumOfKnots-1, 1);
	for(int i=0; i<iNumOfKnots-1;i++)
	{
		gsl_matrix_set (vNormal,i,0,gsl_ran_gaussian(gsl_random_num,1));

		gsl_matrix_set (vBeta,i,0,0);
		for(int j=0; j<=i; j++)
		{
			gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(vNormal,j,0)*gsl_matrix_get(mVar,i,j));
		}

		gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(mMean,i,0));

		sParam.vBeta[i]=gsl_matrix_get(vBeta,i,0);
	}

	/* Seasonal */

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(vBeta,j,0)*gsl_matrix_get(mWTilde,i,j);
		}

	}



	/*print*/
//	for(int j=0; j<iNumOfKnots-1; j++)
//	{
//		cout<<"mean temp at "<<j<<" is "<<gsl_matrix_get(mMeanTemp,j,0)<<endl;
//	}
//	for(int j=0; j<iNumOfKnots-1; j++)
//	{
//		cout<<"mean at "<<j<<" is "<<gsl_matrix_get(mMean,j,0)<<endl;
//	}
//
//
//	for(int j=0; j<iNumOfKnots-1; j++)
//	{
//		cout<<"beta at "<<j<<" is "<<gsl_matrix_get(vBeta,j,0)<<endl;
//	}
//
//
//
//	for(int i=4764; i<4768;  i++)
//	{
//
//		printf ("vSesonal(%d) = %g\n", i, sParam.vS[i]);
//
//	}




	gsl_matrix_free (mYBar);
	gsl_matrix_free (mWBar);
	gsl_matrix_free (mWW);
	gsl_matrix_free (mWY);
	gsl_matrix_free (mMean);
	gsl_matrix_free (mMeanTemp);
	gsl_matrix_free (mVar);
	gsl_matrix_free (mVarInv);

	delete dS1;
	delete dS2;
	delete dY2;
	delete dY1;
}


void DrawSesonalPoly(struct AllParam sParam,int iNumOfObs, double * vTimes,  const gsl_rng * gsl_random_num)
{
	int iNumOfKnots=3;
	unsigned int iTimeInSec=23400;
//	double dDuration=(double)iTimeInSec/iNumOfObs;

	int iNumOfOLSObs=0;
	for(int i=0; i<iNumOfObs;i++)
	{
		if(sParam.vN[i]==0)
		{
			iNumOfOLSObs=iNumOfOLSObs+1;
		}
		else
		{
			iNumOfOLSObs=iNumOfOLSObs+2;
		}
	}
//	cout << "iNumOfOLSObs "<<iNumOfOLSObs<<endl;
	gsl_matrix * mWTilde=gsl_matrix_alloc (iNumOfOLSObs, 2);
	gsl_matrix * mYBar=gsl_matrix_alloc (iNumOfOLSObs, 1);
	gsl_matrix * mWBar=gsl_matrix_alloc(iNumOfOLSObs,iNumOfKnots-1);

	int iTempCount=0;
	for(int i= 0 ; i<iNumOfObs; i++)
	{
//		cout << "vTimes[i] "<<vTimes[i]<<endl;
//		cout << "sParam.vN[i] "<<sParam.vN[i]<<endl;
		if(sParam.vN[i]==0)
		{
			gsl_matrix_set(mWTilde,iTempCount,0, pow(vTimes[i],2)-pow( (double) iTimeInSec, 2)/3 );
			gsl_matrix_set(mWTilde,iTempCount,1,vTimes[i]- ((double) iTimeInSec )/2 );

			iTempCount=iTempCount+1;

		}
		else
		{
			gsl_matrix_set(mWTilde,iTempCount,0, pow(vTimes[i],2)-pow( (double) iTimeInSec, 2)/3 );
			gsl_matrix_set(mWTilde,iTempCount,1,vTimes[i]- ((double) iTimeInSec )/2 );
			gsl_matrix_set(mWTilde,iTempCount+1,0, pow(vTimes[i],2)-pow( (double) iTimeInSec, 2)/3 );
			gsl_matrix_set(mWTilde,iTempCount+1,1,vTimes[i]- ((double) iTimeInSec )/2 );

			iTempCount=iTempCount+2;
		}
//		cout << "iTempCount "<<iTempCount<<endl;
	}

//	for(int i=0; i<iNumOfOLSObs;i++)
//	{
//		for(int j=0; j<iNumOfKnots-1;j++)
//		{
//			cout << "WTilde "<< i<<" "<<j <<" "<<gsl_matrix_get(mWTilde,i,j)<<endl;
//		}
//	}




	double * dY1=new double;
	double * dY2=new double;

	double * dS1=new double;
	double * dS2=new double;

	iTempCount=0;
	for(int i=0; i<iNumOfObs; i++)
	{




		if(sParam.vN[i]==0)
		{
			dY1[0]=sParam.mAuxY[i*2]-(sParam.dMu[0]+sParam.vX[i]);
			dS1[0]=pow(sParam.mAuxH[i*4],0.5);

//			cout<<"XXXXXXXXXXXXXXX"<<endl;
//			cout<<"dY1[0] "<<dY1[0]<<endl;
//			cout<<"dS1[0] "<<dS1[0]<<endl;

			gsl_matrix_set (mYBar,iTempCount,0, dY1[0]/dS1[0] );
			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,iTempCount,j, gsl_matrix_get (mWTilde, iTempCount, j)/dS1[0]);
			}

			iTempCount=iTempCount+1;
		}
		else
		{
			dY1[0]=sParam.mAuxY[i*2]-(sParam.dMu[0]+sParam.vX[i]);
			dY2[0]=sParam.mAuxY[i*2+1]-(sParam.dMu[0]+sParam.vX[i]);
			dS1[0]=pow(sParam.mAuxH[i*4],0.5);
			dS2[0]=pow(sParam.mAuxH[i*4+3],0.5);

//			cout<<"XXXXXXXXXXXXXXX"<<endl;
//			cout<<"dY1[0] "<<dY1[0]<<endl;
//			cout<<"dS1[0] "<<dS1[0]<<endl;
//			cout<<"dY2[0] "<<dY2[0]<<endl;
//			cout<<"dS2[0] "<<dS2[0]<<endl;


			gsl_matrix_set (mYBar,iTempCount,0, dY1[0]/dS1[0]);
			gsl_matrix_set (mYBar,iTempCount+1,0, dY2[0]/dS2[0]);

			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,iTempCount,j, gsl_matrix_get (mWTilde, iTempCount, j)/dS1[0]);
				gsl_matrix_set (mWBar,iTempCount+1,j, gsl_matrix_get (mWTilde, iTempCount+1, j)/dS2[0]);
			}

			iTempCount=iTempCount+2;
		}



	}

//
//	for(int i=0; i<iNumOfOLSObs;i++)
//	{
//		for(int j=0; j<iNumOfKnots-1;j++)
//		{
//			cout << "WBar "<< i<<" "<<j <<" "<<gsl_matrix_get(mWBar,i,j)<<endl;
//		}
//	}
//
//
//	for(int i=0; i<iNumOfOLSObs;i++)
//	{
//
//			cout << "YBar "<< i << " "<<gsl_matrix_get(mYBar,i,0)<<endl;
//
//	}

	/* Variance */
	gsl_matrix * mWW=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVar=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVarInv=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mWBar,0.0, mWW);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			if(j==i)
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j)+sParam.dPriorSeasonalVar[0]);
			}
			else
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j));
			}
		}
	}

	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots-1);
	int s;
	gsl_linalg_LU_decomp (mVarInv, mPermutation, &s);
	gsl_linalg_LU_invert (mVarInv, mPermutation, mVar);



	/* Mean */
	gsl_matrix * mWY=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMean=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMeanTemp=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mYBar,0.0, mWY);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		gsl_matrix_set (mMeanTemp,i,0,gsl_matrix_get(mWY,i,0)+sParam.dPriorSeasonalMean[0]/sParam.dPriorSeasonalVar[0]);
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mVar,mMeanTemp,0.0, mMean);


	/* Draw Normal */
	gsl_linalg_cholesky_decomp (mVar);
	gsl_matrix * vNormal=gsl_matrix_alloc (iNumOfKnots-1, 1);


	/* Draw Beta */
	gsl_matrix * vBeta=gsl_matrix_alloc (iNumOfKnots-1, 1);
	for(int i=0; i<iNumOfKnots-1;i++)
	{
		gsl_matrix_set (vNormal,i,0,gsl_ran_gaussian(gsl_random_num,1));

		gsl_matrix_set (vBeta,i,0,0);
		for(int j=0; j<=i; j++)
		{
			gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(vNormal,j,0)*gsl_matrix_get(mVar,i,j));
		}

		gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(mMean,i,0));
	}
	/* Seasonal */
	iTempCount=0;
	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(vBeta,j,0)*gsl_matrix_get(mWTilde,iTempCount,j);
		}

		if(sParam.vN[i]==0)
		{
			iTempCount=iTempCount+1;
		}
		else
		{
			iTempCount=iTempCount+2;
		}

	}

	gsl_matrix_free (mYBar);
	gsl_matrix_free (mWBar);
	gsl_matrix_free (mWW);
	gsl_matrix_free (mWY);
	gsl_matrix_free (mMean);
	gsl_matrix_free (mMeanTemp);
	gsl_matrix_free (mVar);
	gsl_matrix_free (mVarInv);

	delete dY2;
	delete dY1;
	delete dS2;
	delete dS1;
}


void DrawSesonalTest2(struct AllParam sParam,  int iNumOfObs,const gsl_rng * gsl_random_num)
{
	for(int i=0; i<iNumOfObs; i++)
	{

		if(i % 2 ==0)
		{
			cout <<"even "<<endl;
			sParam.mAuxY[2*i]=sParam.vS[i]+(sParam.dMu[0]+sParam.vX[i])+gsl_ran_gaussian(gsl_random_num,sqrt(0.5) );
			sParam.mAuxY[2*i+1]=sParam.vS[i]+(sParam.dMu[0]+sParam.vX[i])+gsl_ran_gaussian(gsl_random_num,sqrt(2) );

			sParam.mAuxH[4*i]=0.5;
			sParam.mAuxH[4*i+1]=0;
			sParam.mAuxH[4*i+2]=0;
			sParam.mAuxH[4*i+3]=2;
		}
		else
		{
			cout <<"odd "<<endl;
			sParam.mAuxY[2*i]=sParam.vS[i]+(sParam.dMu[0]+sParam.vX[i])+gsl_ran_gaussian(gsl_random_num,sqrt(3) );
			sParam.mAuxY[2*i+1]=sParam.vS[i]+(sParam.dMu[0]+sParam.vX[i])+gsl_ran_gaussian(gsl_random_num,sqrt(0.3) );

			sParam.mAuxH[4*i]=3;
			sParam.mAuxH[4*i+1]=0;
			sParam.mAuxH[4*i+2]=0;
			sParam.mAuxH[4*i+3]=0.3;
		}
	}
}

void DrawSesonalTest(struct AllParam sParam,double *vYBar, int iNumOfObs, int iNumOfKnots,gsl_matrix * mWBar ,const gsl_rng * gsl_random_num)
{
	gsl_matrix * mYBar=gsl_matrix_alloc (iNumOfObs, 1);


	for(int i=0; i<iNumOfObs; i++)
	{

		gsl_matrix_set (mYBar,i,0,vYBar[i] );

	}

	/* Variance */
	gsl_matrix * mWW=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVar=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVarInv=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mWBar,0.0, mWW);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			if(j==i)
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j)+sParam.dPriorSeasonalVar[0]);
			}
			else
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j));
			}
		}
	}

	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots-1);
	int s;
	gsl_linalg_LU_decomp (mVarInv, mPermutation, &s);
	gsl_linalg_LU_invert (mVarInv, mPermutation, mVar);


	/* Mean */
	gsl_matrix * mWY=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMean=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMeanTemp=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mYBar,0.0, mWY);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		gsl_matrix_set (mMeanTemp,i,0,gsl_matrix_get(mWY,i,0)+sParam.dPriorSeasonalMean[0]/sParam.dPriorSeasonalVar[0]);
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mVar,mMeanTemp,0.0, mMean);


//		for(int i=0; i<iNumOfKnots-1; i++)
//		{
//			for(int j=0; j<iNumOfKnots-1; j++)
//			{
//				 printf ("mWW(%d,%d) = %g\n", i, j, gsl_matrix_get (mWW, i, j));
//			}
//		}
//
//		for(int i=0; i<iNumOfKnots-1; i++)
//		{
//			printf ("mWY(%d,%d) = %g\n", i, 0, gsl_matrix_get (mWY, i,0));
//		}
//
//
//		for(int i=0; i<iNumOfKnots-1; i++)
//				{
//					for(int j=0; j<iNumOfKnots-1; j++)
//					{
//						 printf ("mVar(%d,%d) = %g\n", i, j, gsl_matrix_get (mVar, i, j));
//					}
//				}
//
//				for(int i=0; i<iNumOfKnots-1; i++)
//				{
//					printf ("mMean(%d,%d) = %g\n", i, 0, gsl_matrix_get (mMean, i,0));
//				}
	/* Draw Normal */
	gsl_linalg_cholesky_decomp (mVar);
	gsl_matrix * vNormal=gsl_matrix_alloc (iNumOfKnots-1, 1);


//	for(int i=0; i<iNumOfKnots-1; i++)
//				{
//					for(int j=0; j<iNumOfKnots-1; j++)
//					{
//						 printf ("mChol(%d,%d) = %g\n", i, j, gsl_matrix_get (mVar, i, j));
//					}
//				}

	/* Draw Beta */
	gsl_matrix * vBeta=gsl_matrix_alloc (iNumOfKnots-1, 1);
	for(int i=0; i<iNumOfKnots-1;i++)
	{
		gsl_matrix_set (vNormal,i,0,gsl_ran_gaussian(gsl_random_num,1));
		cout<<"Normal "<<i<< " is "<<gsl_matrix_get(vNormal,i,0) <<endl;
		gsl_matrix_set (vBeta,i,0,0);
		for(int j=0; j<=i; j++)
		{
			gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(vNormal,j,0)*gsl_matrix_get(mVar,i,j));
		}

		gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(mMean,i,0));

	}

	/* Seasonal */
	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(vBeta,j,0)*gsl_matrix_get(mWBar,i,j);
		}
	}

	gsl_matrix_free (mYBar);
	gsl_matrix_free (mWW);
	gsl_matrix_free (mWY);
	gsl_matrix_free (mMean);
	gsl_matrix_free (mMeanTemp);
	gsl_matrix_free (mVar);
	gsl_matrix_free (mVarInv);

}

int Cholesky(double* mA,  int iSize, double* mL)
{
	gsl_matrix * mM = gsl_matrix_alloc (iSize, iSize);
	for(int i=0;i<iSize;i++)
	{
		for(int j=0;j<iSize;j++)
		{
			gsl_matrix_set (mM, i, j, mA[i*iSize+j]);
		}
	}

	gsl_set_error_handler_off();
	int status =gsl_linalg_cholesky_decomp (mM);
	gsl_set_error_handler (NULL);
	int iFlag;
	iFlag=1;
	for(int i=0;i<iSize;i++)
	{
		for(int j=0;j<iSize;j++)
		{
			if(i<=j)
			{
				mL[i*iSize+j]=gsl_matrix_get (mM, i, j);
			}
			else
			{
				mL[i*iSize+j]=0;
			}

			if(i==j)
			{
				if((mL[i*iSize+j]!=mL[i*iSize+j]) | (mL[i*iSize+j]<0))
				{
					iFlag=0;
				}
			}
		}
	}

	if(status)
	{
		 iFlag=0;
	}
	gsl_matrix_free (mM );

	return iFlag;
}

void CaculateDNBLL(const double * vData,  int iNumOfObs,  struct AllParam sParam,double * vDNBLL )
{
	double* dLambda=new double;
	double* dLambdaRatio=new double;
	double* dNuRatio=new double;
	int * iAbsX=new int;
	double * dDNB=new double;
	vDNBLL[0]=0;


	for(int i=0; i<iNumOfObs; i++)
	{

		dLambda[0]=exp(sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]);
		dLambdaRatio[0]=dLambda[0]/(dLambda[0]+sParam.dNu[0]);
		dNuRatio[0]=1-dLambdaRatio[0];
		iAbsX[0]=abs(vData[i]);

		dDNB[0]=exp(2*sParam.dNu[0]*log(dNuRatio[0])+iAbsX[0]*log(dLambdaRatio[0])+lgamma(sParam.dNu[0]+iAbsX[0])-lgamma(sParam.dNu[0])-lgamma(iAbsX[0]+1)
					+log(gsl_sf_hyperg_2F1(sParam.dNu[0]+iAbsX[0],sParam.dNu[0], iAbsX[0]+1 , pow(dLambdaRatio[0],2))));

		if(vData[i]==0)
		{
			vDNBLL[0]=vDNBLL[0]+log(sParam.dGamma[0]+(1-sParam.dGamma[0])*dDNB[0]) ;
		}
		else
		{
			vDNBLL[0]=vDNBLL[0]+log((1-sParam.dGamma[0])*dDNB[0]) ;
		}


	}

	delete dLambda;
	delete dLambdaRatio;
	delete dNuRatio;
	delete iAbsX;
	delete dDNB;
}

void NuInBlockLogPosteriorDNB(double * dLogNu, double *vData, struct AllParam sParam,  double *dPriorA, double *dPriorB, unsigned int iNumOfObs,
	double * dLogPosterior)
{

	double dNu=exp(dLogNu[0]);
	double dOldNu=sParam.dNu[0];
	sParam.dNu[0]=dNu;



	CaculateDNBLL(vData, iNumOfObs, sParam,dLogPosterior );

	// dLogPosterior[0]+=dPriorA[0]*log(dPriorB[0]) +dPriorA[0]*log(dNu) -dPriorB[0]*dNu-lgamma(dPriorA[0]);
	dLogPosterior[0]+=dPriorA[0]*log(dPriorB[0]) + (dPriorA[0]-1)*log(dNu) -dPriorB[0]*dNu-lgamma(dPriorA[0]);
	sParam.dNu[0]=dOldNu;

}

void DrawNuInBlockAdaptiveRWDNB( int iNumOfObs,  int iNumOfIter, double * vData, const gsl_rng * gsl_random_num, struct AllParam sParam, double* dCovar,  double* dSum  )
{


	/* Proposal */
	double dOmega1=0.05;
	double dU=gsl_rng_uniform(gsl_random_num);
	double dNewNuLog;
	double dSigma;
	double * dLogPosteriorNew=new double;
	double * dLogPosteriorOld=new double;

	if(iNumOfIter<=1)
	{
		dSigma=0.1;
	}
	else if(iNumOfIter==2)
	{
		dSigma=pow(dCovar[0]-0.5*dSum[0]*dSum[0],0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}
	else
	{
		dSigma=pow(dCovar[0], 0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	}

	if(dU<=dOmega1)
	{
		/* Static proposal */
//		dNewNuLog=log(sParam.dNu[0])+gsl_ran_gaussian(gsl_random_num,0.1);
		dNewNuLog=log(sParam.dNu[0])+gsl_ran_tdist(gsl_random_num,5)*0.1;
	}
	else
	{
		/* Proposal based on empirical covariance */
//		dNewNuLog=log(sParam.dNu[0])+gsl_ran_gaussian(gsl_random_num,dSigma*2.38);
		dNewNuLog=log(sParam.dNu[0])+gsl_ran_tdist(gsl_random_num,5)*dSigma*2.38;
	}


	/* Acceptance Rate */

	double dOldNuLog=log(sParam.dNu[0]);


	NuInBlockLogPosteriorDNB(&dNewNuLog, vData, sParam,  sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorNew);
	NuInBlockLogPosteriorDNB(&dOldNuLog, vData, sParam,  sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorOld);


//	cout <<"dLogPosteriorNew[0] "<<dLogPosteriorNew[0]<<endl;
//	cout <<"dLogPosteriorOld[0] "<<dLogPosteriorOld[0]<<endl;

	double dLogU=log(gsl_rng_uniform(gsl_random_num));
    double dAccept=fmin(dLogPosteriorNew[0]-dLogPosteriorOld[0],0);

    cout <<"dAccpet "<<dAccept<<" dLogU "<< dLogU<<endl;
	if(dLogU<=dAccept)
	{
		sParam.dNu[0]=exp(dNewNuLog);
//		cout <<"Accepted! "<<endl;
	}



	/* Covariance */
	if(iNumOfIter<=1)
	{
		dCovar[0]=dCovar[0]+pow(log(sParam.dNu[0]),2);
		dSum[0]=dSum[0]+log(sParam.dNu[0]);

	}
	else if(iNumOfIter==2)
	{
		dCovar[0]=(iNumOfIter-1)*pow(dSigma,2)/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(log(sParam.dNu[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+log(sParam.dNu[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));


	}
	else
	{
		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(log(sParam.dNu[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+log(sParam.dNu[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));
	}

//	cout << "dSum " << dSum[0]<< endl;
//	cout << "dCovar " << dCovar[0]<< endl;
//	cout << "proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
//	cout << "accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;

}

void DrawNuInBlockDiscreteDNB(double * vData, struct AllParam sParam ,int iNumOfObs, const gsl_rng * gsl_random_num, double dResolution,  double dDelta)
{
	/* Propose draw */
	double dU=gsl_rng_uniform(gsl_random_num);
	int iN=2*dDelta/dResolution+1;
	double * dLogPosteriorNew=new double;
	double * dLogPosteriorOld=new double;
	double dCum=1.0/iN;
	double dNuProposed=sParam.dNu[0]-dDelta;
	while(dCum<dU)
	{
		dCum=dCum+1.0/iN;
		dNuProposed=dNuProposed+dResolution;
	}

	cout<< "dNuProposed "<<dNuProposed<<endl;
	cout<< "dNuCurrent "<<sParam.dNu[0]<<endl;

	double dNewNuLog=log(dNuProposed);
	double dOldNuLog=log(sParam.dNu[0]);

	// if((dNuProposed <= 128) & (dNuProposed >=0.1))
	if((dNuProposed <= 128.0) & (dNuProposed >=2.0))
	{
		/* Calculate acceptance probability */


		NuInBlockLogPosteriorDNB(&dNewNuLog, vData, sParam,  sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorNew);
		NuInBlockLogPosteriorDNB(&dOldNuLog, vData, sParam,  sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorOld);


		double dAccept=fmin(dLogPosteriorNew[0]- dLogPosteriorOld[0],0);

		/* Accept or reject*/
		double dLogU=log(gsl_rng_uniform(gsl_random_num));

		cout << "dAccept "<<dAccept<<" dLogU "<< dLogU<<endl;
		if(dLogU<=dAccept)
		{
			sParam.dNu[0]=dNuProposed;
		}
	}

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;

}

void SystemMatricesWithSeasonal( int iNumOfKnots, struct AllParam sParam,   double* mZ, double* mT, double *mQ, double* mA1, double* mP1)
{
	/* Z */
	for(int i=0; i<iNumOfKnots-1;i++)
	{
		mZ[i]=0;
		mZ[(iNumOfKnots-1+2)+i]=0;
	}
	for(int i=0; i<2;i++)
	{
		mZ[iNumOfKnots-1+i]=1;
		mZ[(iNumOfKnots-1+2)+iNumOfKnots-1+i]=1;
	}





	/* T */
	for(int i=0;i<(iNumOfKnots-1+2);i++)
	{
		for(int j=0;j<(iNumOfKnots-1+2);j++)
		{
			if( (i==j) & (i<iNumOfKnots-1+1))
			{
				mT[i*(iNumOfKnots-1+2)+j ]=1;
			}
			else if ( (i==j) & (i==iNumOfKnots-1+1))
			{
				mT[i*(iNumOfKnots-1+2)+j ]=sParam.dPhi[0];
			}
			else
			{
				mT[i*(iNumOfKnots-1+2)+j ]=0;
			}
		}
	}




	/* Q */
	for(int i=0;i<(iNumOfKnots-1+2);i++)
	{
		for(int j=0;j<(iNumOfKnots-1+2);j++)
		{
			if ( (i==j) & (i==iNumOfKnots-1+1))
			{
				mQ[i*(iNumOfKnots-1+2)+j ]=sParam.dSigma2[0];
			}
			else
			{
				mQ[i*(iNumOfKnots-1+2)+j ]=0;
			}
		}
	}




	/* A1 */
	for(int i=0;i<(iNumOfKnots-1);i++)
	{
		mA1[i]=sParam.dPriorSeasonalMean[0];
	}

	mA1[iNumOfKnots-1]=sParam.dPriorMuMean[0];
	mA1[iNumOfKnots]=0;



	/* P1 */
	for(int i=0;i<(iNumOfKnots-1+2);i++)
	{
		for(int j=0;j<(iNumOfKnots-1+2);j++)
		{
			if( (i==j) & (i<iNumOfKnots-1))
			{
				mP1[i*(iNumOfKnots-1+2)+j ]=sParam.dPriorSeasonalVar[0];
			}
			else if ( (i==j) & (i==iNumOfKnots-1))
			{
				mP1[i*(iNumOfKnots-1+2)+j ]=sParam.dPriorMuSigma2[0];
			}
			else if ( (i==j) & (i==iNumOfKnots-1+1))
			{
				mP1[i*(iNumOfKnots-1+2)+j ]=sParam.dSigma2[0]/(1-sParam.dPhi[0]*sParam.dPhi[0]);
			}
			else
			{
				mP1[i*(iNumOfKnots-1+2)+j ]=0;
			}
		}
	}


//	for(int i=0; i<(iNumOfKnots-1+2)*2; i++)
//	{
//		cout<<"mZ at "<<i<<" is "<<mZ[i]<<endl;
//	}
//	for(int i=0; i<(iNumOfKnots-1+2)*(iNumOfKnots-1+2); i++)
//	{
//		cout<<"mT at "<<i<<" is "<<mT[i]<<endl;
//	}
//
//	for(int i=0; i<(iNumOfKnots-1+2)*(iNumOfKnots-1+2); i++)
//	{
//		cout<<"mQ at "<<i<<" is "<<mQ[i]<<endl;
//	}
//
//	for(int i=0; i<(iNumOfKnots-1+2); i++)
//	{
//		cout<<"mA1 at "<<i<<" is "<<mA1[i]<<endl;
//	}
//
//	for(int i=0; i<(iNumOfKnots-1+2)*(iNumOfKnots-1+2); i++)
//	{
//		cout<<"mP1 at "<<i<<" is "<<mP1[i]<<endl;
//	}




}

void KalmanFilterWithSeasonal(int iNumOfKnots,   gsl_matrix * mWTilde , double* mY, int* mN, double* mH, double* mZ, double* mT, double *mQ,
		double* mA1, double* mP1, int iDimOfObs, int iDimOfStates,  int iNumOfObs, double* mA, double * mP, double* mFInv, double* mV, double* mL)
{

	double * mZTrans= new double[iDimOfObs*iDimOfStates];


	for(int d=0; d<iDimOfStates; d++)
	{
		mA[d]=mA1[d];
		for(int k=0; k<iDimOfStates; k++)
		{
			mP[d*iDimOfStates+k]=mP1[d*iDimOfStates+k];
//			cout<<"mP[0] "<< d*iDimOfStates+k <<" " <<mP[d*iDimOfStates+k]<<endl;
		}
	}

	double* mNextA;
	double* mNextP;

	double * mZa= new double[iDimOfObs];
	double * mZPZ= new double[iDimOfObs*iDimOfObs];
	double * mK= new double[iDimOfObs*iDimOfStates];
	double * mPZFInv = new double[iDimOfStates*iDimOfObs];
	double * mZFInv = new double[iDimOfStates*iDimOfObs];
	double * mF= new double[iDimOfObs*iDimOfObs];
	double * mKZ= new double[iDimOfStates*iDimOfStates];
	double * mTa= new double[iDimOfStates];
	double * mKv= new double[iDimOfStates];
	double * mLTrans=new double[iDimOfStates*iDimOfStates];
	double * mPL=new double[iDimOfStates*iDimOfStates];
	double * mTPL=new double[iDimOfStates*iDimOfStates];

	for(int i=0; i<iNumOfObs; i++)
	{
		/* Seting the time varying Z matrix */
		for(int k=0; k<iNumOfKnots-1;k++)
		{
				mZ[k]=gsl_matrix_get(mWTilde,i,k);
				mZ[(iNumOfKnots-1+2)+k]=gsl_matrix_get(mWTilde,i,k);
		}
//		for(int l=0; l<(iNumOfKnots-1+2)*2; l++)
//		{
//				cout<< "At obs "<<i <<" mZ at "<<l<<" is "<<mZ[l]<<endl;
//		}

		/* Setting the observation dimension to 1 when there is only one observation */
		if(mN[i]==0)
		{
			iDimOfObs=1;
		}

		MatrixTrans( mZ, iDimOfObs, iDimOfStates, mZTrans);

//		for(int l=0; l<iDimOfStates*iDimOfObs; l++)
//		{
//			cout<< "At obs "<<i <<" mZ at "<<l<<" is "<<mZTrans[l]<<endl;
//		}


		/* v_t */
		if(i==0)
		{
			MatrixMulti(mZ, mA1,iDimOfObs,iDimOfStates, iDimOfStates, 1, mZa);
		}
		else
		{
			MatrixMulti(mZ, mA,iDimOfObs,iDimOfStates, iDimOfStates, 1, mZa);
		}

		MatrixSub(mY, mZa, iDimOfObs, 1, mV);



		/* F_t */
		if(i==0)
		{
			MatrixSandwitch(mZTrans, mP1,iDimOfStates,iDimOfObs, iDimOfStates, iDimOfStates, mZPZ);
		}
		else
		{
			MatrixSandwitch(mZTrans, mP,iDimOfStates,iDimOfObs, iDimOfStates, iDimOfStates, mZPZ);
		}

		MatrixAdd(mZPZ, mH, iDimOfObs, iDimOfObs, mF);

//		for(int k=0; k<iDimOfObs; k++)
//		{
//			cout <<"At obs "<<i << " F  "<<k<<" is "<<mF[iDimOfObs*k+k]<<endl;
//		}

		/* K_t */
		MatirxInv2x2(mF , iDimOfObs,mFInv);
		MatrixMulti(mZTrans, mFInv,iDimOfStates, iDimOfObs, iDimOfObs, iDimOfObs, mZFInv);
		MatrixMulti(mP, mZFInv,iDimOfStates, iDimOfStates, iDimOfStates, iDimOfObs, mPZFInv);
		MatrixMulti(mT, mPZFInv,iDimOfStates, iDimOfStates, iDimOfStates, iDimOfObs, mK);

//		for(int k=0; k<iDimOfObs*iDimOfStates; k++)
//		{
//			cout <<"At obs "<<i << " K  "<<k<<" is "<<mK[k]<<endl;
//		}

		/* L_t */
		MatrixMulti(mK, mZ,iDimOfStates, iDimOfObs, iDimOfObs, iDimOfStates, mKZ);
//		for(int k=0; k<iDimOfStates*iDimOfStates; k++)
//		{
//			cout <<"At obs "<<i << " KZ  "<<k<<" is "<<mKZ[k]<<endl;
//		}

		MatrixSub(mT, mKZ, iDimOfStates, iDimOfStates, mL);

//		for(int k=0; k<iDimOfStates*iDimOfStates; k++)
//		{
//			cout <<"At obs "<<i << " L  "<<k<<" is "<<mL[k]<<endl;
//		}

		if(i<iNumOfObs-1)
		{
			mNextA=mA+iDimOfStates;
			mNextP=mP+iDimOfStates*iDimOfStates;
			/* A_t */
			MatrixMulti(mT, mA,iDimOfStates, iDimOfStates, iDimOfStates, 1, mTa);
//			for(int k=0; k<iDimOfStates; k++)
//			{
//				cout <<"At obs "<<i << " Ta  "<<k<<" is "<<mTa[k]<<endl;
//			}
			MatrixMulti(mK, mV,iDimOfStates, iDimOfObs, iDimOfObs, 1, mKv);
//			for(int k=0; k<iDimOfStates; k++)
//			{
//							cout <<"At obs "<<i << " mKv  "<<k<<" is "<<mKv[k]<<endl;
//			}

			MatrixAdd(mTa, mKv, iDimOfStates, 1, mNextA);
//			for(int k=0; k<iDimOfStates; k++)
//			{
//				cout <<"At obs "<<i << " mNextA  "<<k<<" is "<<mNextA[k]<<endl;
//			}


			/* P_t */
			MatrixTrans( mL, iDimOfStates, iDimOfStates, mLTrans);
			MatrixMulti(mP, mLTrans,iDimOfStates, iDimOfStates, iDimOfStates,  iDimOfStates, mPL);
			MatrixMulti(mT, mPL,iDimOfStates, iDimOfStates, iDimOfStates,  iDimOfStates, mTPL);
			MatrixAdd(mTPL, mQ, iDimOfStates, iDimOfStates, mNextP);

//			cout<<"mA[0] "<<mA[0]<<endl;
//			cout<<"mA[1] "<<mA[1]<<endl;


			/* Set Dim back */
			if(mN[i]==0)
			{
				iDimOfObs=2;
			}

			/* Update pointers */
			mA=mA+iDimOfStates;
			mP=mP+iDimOfStates*iDimOfStates;
			mV=mV+iDimOfObs;
			mFInv=mFInv+iDimOfObs*iDimOfObs;
			mL=mL+iDimOfStates*iDimOfStates;
			mH=mH+iDimOfObs*iDimOfObs;
			mY=mY+iDimOfObs;
		}

//		cout<<"XXXXXXXXXXXXXXXXX"<<endl;

	}

	delete [] mZTrans;
	delete [] mZa;
	delete [] mZPZ;
	delete [] mK;
	delete [] mPZFInv;
	delete [] mZFInv;
	delete [] mF;
	delete [] mTa;
	delete [] mKv;
	delete [] mPL;
	delete [] mTPL;
	delete [] mKZ;
	delete [] mLTrans;
}

void KalmanSmootherWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,int * vN, double* mZ,  double * mA, double *mP, double *mV,
		double* mFInv, double* mL,   int iDimOfObs,int iDimOfStates,  int iNumOfObs, double *mAHat, double* mVHat )
{
	double * mZTrans= new double[iDimOfObs*iDimOfStates];

	double * mR= new double[iDimOfStates];
	double * mNextR= new double[iDimOfStates];
	double * mN= new double[iDimOfStates*iDimOfStates];
	double * mNextN = new double[iDimOfStates*iDimOfStates];

	double * mLTrans=new double[iDimOfStates*iDimOfStates];
	double * mFv =new double[iDimOfObs];
	double * mZFv=new double[iDimOfStates];
	double * mLR=new double[iDimOfStates];
	double * mPR=new double[iDimOfStates];
	double * mZFZ= new double[iDimOfStates*iDimOfStates];
	double * mLNL= new double[iDimOfStates*iDimOfStates];
	double * mPNP= new double[iDimOfStates*iDimOfStates];

	for(int i=0; i<iDimOfStates; i++)
	{
		mR[i]=0;
		for(int j=0; j<iDimOfStates; j++)
		{
			mN[i*iDimOfStates+j]=0;
		}
	}

	/* Set the pointer to the end of the arrays (backward recursion)*/
	mP=mP+(iNumOfObs-1)*iDimOfStates*iDimOfStates;
	mAHat=mAHat+(iNumOfObs-1)*iDimOfStates;
	mVHat=mVHat+(iNumOfObs-1)*iDimOfStates*iDimOfStates;
	mA=mA+(iNumOfObs-1)*iDimOfStates;
	mV=mV+(iNumOfObs-1)*iDimOfObs;
	mFInv=mFInv+(iNumOfObs-1)*iDimOfObs*iDimOfObs;
	mL=mL+(iNumOfObs-1)*iDimOfStates*iDimOfStates;

	for(int i=iNumOfObs-1; i>-1; i--)
	{

		for(int k=0; k<iNumOfKnots-1;k++)
		{
						mZ[k]=gsl_matrix_get(mWTilde,i,k);
						mZ[(iNumOfKnots-1+2)+k]=gsl_matrix_get(mWTilde,i,k);
		}

		if(vN[i]==0)
		{
			iDimOfObs=1;
		}

		MatrixTrans( mZ, iDimOfObs, iDimOfStates, mZTrans);
		/* mR */
		MatrixTrans( mL, iDimOfStates, iDimOfStates, mLTrans);
		MatrixMulti(mFInv, mV,iDimOfObs, iDimOfObs, iDimOfObs, 1, mFv);
		MatrixMulti(mZTrans, mFv,iDimOfStates, iDimOfObs, iDimOfObs, 1, mZFv);
		MatrixMulti(mLTrans, mR,iDimOfStates, iDimOfStates, iDimOfStates,1, mLR);
		MatrixAdd( mZFv, mLR, iDimOfStates, 1, mNextR);

		/* AHat */
		MatrixMulti(mP, mNextR,iDimOfStates, iDimOfStates, iDimOfStates,1, mPR);
		MatrixAdd( mA, mPR, iDimOfStates, 1, mAHat);

		/* mN*/
		MatrixSandwitch(mZ, mFInv,iDimOfObs,iDimOfStates, iDimOfObs, iDimOfObs, mZFZ);
		MatrixSandwitch(mL, mN,iDimOfStates,iDimOfStates, iDimOfStates, iDimOfStates, mLNL);
		MatrixAdd( mZFZ, mLNL, iDimOfStates, iDimOfStates, mNextN);

		/* VHat */
		MatrixSandwitch(mP, mNextN,iDimOfStates,iDimOfStates, iDimOfStates, iDimOfStates, mPNP);
		MatrixSub( mP, mPNP, iDimOfStates, iDimOfStates, mVHat);


		if(vN[i]==0)
		{
			iDimOfObs=2;
		}

		/* Update R and N values */
		for(int k=0; k<iDimOfStates; k++)
		{
			mR[k]=mNextR[k];
			for(int j=0; j<iDimOfStates; j++)
			{
				mN[k*iDimOfStates+j]=mNextN[k*iDimOfStates+j];
			}
		}
		if(i>0)
		{
			/* Update pointers */
			mP=mP-iDimOfStates*iDimOfStates;
			mAHat=mAHat-iDimOfStates;
			mVHat=mVHat-iDimOfStates*iDimOfStates;
			mA=mA-iDimOfStates;
			mV=mV-iDimOfObs;
			mFInv=mFInv-iDimOfObs*iDimOfObs;
			mL=mL-iDimOfStates*iDimOfStates;
		}
	}

	delete [] mZTrans;
	delete [] mR;
	delete [] mNextR;
	delete [] mN;
	delete [] mNextN ;
	delete [] mLTrans;
	delete [] mFv ;
	delete [] mZFv;
	delete [] mLR;
	delete [] mPR;
	delete [] mZFZ;
	delete [] mLNL;
	delete [] mPNP;
}

void CalculateLLWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam,int iNumOfObs,  int iDimOfObs,
		int iDimOfStates, double * dLL )
{
	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatricesWithSeasonal(iNumOfKnots, sParam, mZ, mT, mQ,  mA1,  mP1);

	double * mA =new double[iNumOfObs*iDimOfStates];
	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
	double * mV =new double[iNumOfObs*iDimOfObs];
	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];
	KalmanFilterWithSeasonal(iNumOfKnots,   mWTilde ,sParam.mAuxY, sParam.vN, sParam.mAuxH,   mZ, mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mY at " << i<< " is "<< sParam.mAuxY[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*4; i++)
//	{
//		cout <<"mH at " << i<< " is "<< sParam.mAuxH[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//			cout <<"mV at " << i<< " is "<< mV[i]<< endl;
//	}

//	string sYFile="mY.csv";
//	WriteOutDoubleArray(sParam.mAuxY, iNumOfObs, 2, sYFile);
//	string sHFile="mH.csv";
//	WriteOutDoubleArray(sParam.mAuxH, iNumOfObs, 4, sHFile);
//	string sNFile="mN.csv";
//	WriteOutIntArray(sParam.vN, iNumOfObs, 1, sNFile);



	dLL[0]=0;
//	cout << "LL at init is " << dLL[0] << endl;
	double * dDetFInv= new double;
	double * dvFv=new double;


	for(int i=0; i<iNumOfObs; i++)
	{
		if(sParam.vN[i]==0)
		{
			iDimOfObs=1;
		}

		Determinant2x2(mFInv,iDimOfObs, dDetFInv);
		MatrixSandwitch(mV, mFInv,iDimOfObs,1, iDimOfObs, iDimOfObs, dvFv);
		dLL[0]=dLL[0]-0.5*iDimOfObs*log(2*PI)+0.5*(log(dDetFInv[0])-dvFv[0]);
//		cout << "LL Cont at " << i << " is " << 0.5*(log(dDetFInv[0])-dvFv[0]) << endl;
//		cout << "LL at " << i << " is " << dLL[0] << endl;

		if(sParam.vN[i]==0)
		{
			iDimOfObs=2;
		}
		if(i<iNumOfObs-1)
		{
			mV=mV+iDimOfObs;
			mFInv=mFInv+iDimOfObs*iDimOfObs;
		}
	}


	delete [] mZ;
	delete [] mT;
	delete [] mQ;
	delete [] mA1;
	delete [] mP1 ;
	delete [] mA;
	delete [] mP ;
	mFInv=mFInv-(iNumOfObs-1)*iDimOfObs*iDimOfObs;
	mV=mV-(iNumOfObs-1)*iDimOfObs;
	delete [] mFInv;
	delete [] mV;
	delete [] mL;
	delete dDetFInv;
	delete dvFv;


}

void SSFRecursionWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam, int * vN, double *mH, double* mZ, double* mT, double *mQ, double* mA1, double* mP1,
		int iNumOfObs,  int iDimOfObs, int iDimOfStates, const gsl_rng * gsl_random_num, double* mYPlus, double* mAPlus)
{



	double * mQChol =new double[iDimOfStates*iDimOfStates];
	double * mHChol =new double[iDimOfObs*iDimOfObs];
	double * mP1Chol =new double[iDimOfStates*iDimOfStates];


	for(int i=0;i<(iNumOfKnots-1+2);i++)
	{
		for(int j=0;j<(iNumOfKnots-1+2);j++)
		{
			if ( (i==j) & (i==iNumOfKnots-1+1))
			{
				mQChol[i*(iNumOfKnots-1+2)+j ]=sqrt(sParam.dSigma2[0]);
			}
			else
			{
				mQChol[i*(iNumOfKnots-1+2)+j ]=0;
			}
		}
	}

	for(int i=0;i<(iNumOfKnots-1+2);i++)
	{
		for(int j=0;j<(iNumOfKnots-1+2);j++)
		{
			if( (i==j) & (i<iNumOfKnots-1))
			{
				mP1Chol[i*(iNumOfKnots-1+2)+j ]=sqrt(sParam.dPriorSeasonalVar[0]);
			}
			else if ( (i==j) & (i==iNumOfKnots-1))
			{
				mP1Chol[i*(iNumOfKnots-1+2)+j ]=sqrt(sParam.dPriorMuSigma2[0]);
			}
			else if ( (i==j) & (i==iNumOfKnots-1+1))
			{
				mP1Chol[i*(iNumOfKnots-1+2)+j ]=sqrt(sParam.dSigma2[0]/(1-sParam.dPhi[0]*sParam.dPhi[0]));
			}
			else
			{
				mP1Chol[i*(iNumOfKnots-1+2)+j ]=0;
			}
		}
	}

//	for(int i=0; i<(iNumOfKnots-1+2)*(iNumOfKnots-1+2); i++)
//	{
//		cout<<"mQChol at "<<i<<" is "<<mQChol[i]<<endl;
//	}
//
//	for(int i=0; i<(iNumOfKnots-1+2)*(iNumOfKnots-1+2); i++)
//	{
//		cout<<"mP1Chol at "<<i<<" is "<<mP1Chol[i]<<endl;
//	}


//	cout << "iDimOfObs "<<iDimOfObs<<endl;
//	cout<< "iDimOfStates "<< iDimOfStates<<endl;

	double * mEps =new double[iDimOfObs];
	double * mEta =new double[iDimOfStates];
	double * mEpsDraw =new double[iDimOfObs];
	double * mEtaDraw =new double[iDimOfStates];
	double * mTa =new double[iDimOfStates];
	double * mNextAPlus =new double[iDimOfStates];
	double * mZa =new double[iDimOfObs];
	double * mNextYPlus =new double[iDimOfObs];

	/* Random values to test SSF */
//	double * mRandom=new double[(iDimOfObs+iDimOfStates)*(iNumOfObs+1)];
//	for(int k=0; k<iDimOfStates;k++)
//	{
//		mRandom[(iDimOfObs+iDimOfStates)*iNumOfObs+k]=0;
//	}
//	for(int k=0; k<2;k++)
//	{
//		mRandom[iDimOfStates+k]=0;
//	}


	for(int i=0; i<iNumOfObs; i++)
	{

		/* Seting the time varying Z matrix */
		for(int k=0; k<iNumOfKnots-1;k++)
		{
			mZ[k]=gsl_matrix_get(mWTilde,i,k);
			mZ[(iNumOfKnots-1+2)+k]=gsl_matrix_get(mWTilde,i,k);
		}

		if(vN[i]==0)
		{
			iDimOfObs=1;

		}
//		cout << "SSF Recuriosn H Chol START " <<endl;
//		cout<<"vN[i] "<<vN[i]<<endl;

		Cholesky(mH, iDimOfObs, mHChol);

//		cout<<"mH[0] "<<mH[0]<<endl;
//		cout<<"mH[3] "<<mH[3]<<endl;
//		cout<<"mHChol[0] "<<mHChol[0]<<endl;
//		cout<<"mHChol[1] "<<mHChol[1]<<endl;
//		cout<<"mHChol[2] "<<mHChol[2]<<endl;
//		cout<<"mHChol[3] "<<mHChol[3]<<endl;

//		cout << "SSF Recuriosn H Chol END " <<endl;
		/*Draw random numbers */
		for (int k=0;k<iDimOfObs;k++)
		{
			mEpsDraw[k]=gsl_ran_gaussian(gsl_random_num,1);

		}
		for (int k=0;k<iDimOfStates;k++)
		{
			mEtaDraw[k]=gsl_ran_gaussian(gsl_random_num,1);

		}


		/* Recursion */
		if(i==0)
		{
			MatrixMulti(mP1Chol, mEtaDraw,iDimOfStates, iDimOfStates, iDimOfStates, 1, mEta);
			//cout << "mEta[0] " <<mEta[0] <<endl;
			MatrixAdd(mA1, mEta, iDimOfStates, 1, mNextAPlus);
			//cout << "mAPlus[0] " <<mNextAPlus[0] <<endl;
		}
		else
		{
//			cout<<"mQChol[0] "<<mQChol[0]<<endl;
//			cout<<"mQChol[1] "<<mQChol[1]<<endl;
//			cout<<"mQChol[2] "<<mQChol[2]<<endl;
//			cout<<"mQChol[3] "<<mQChol[3]<<endl;

			MatrixMulti(mQChol, mEtaDraw,iDimOfStates, iDimOfStates, iDimOfStates, 1, mEta);

//			cout<<"Obs "<<i <<" mEtaDraw[0] "<<mEtaDraw[0]<<endl;
//			cout<<"Obs "<<i <<" mEtaDraw[1] "<<mEtaDraw[1]<<endl;
//			cout<<"Obs "<<i <<" mEta[0] "<<mEta[0]<<endl;
//			cout<<"Obs "<<i <<" mEta[1] "<<mEta[1]<<endl;

			MatrixMulti(mT, mAPlus,iDimOfStates, iDimOfStates, iDimOfStates, 1,mTa);
			MatrixAdd(mTa, mEta, iDimOfStates, 1, mNextAPlus);
		}

//		cout<<iDimOfObs<<endl;
//		cout<<"mH[0] "<<mH[0]<<endl;
//		cout<<"mH[1] "<<mH[1]<<endl;
//		cout<<"mH[2] "<<mH[2]<<endl;
//		cout<<"mH[3] "<<mH[3]<<endl;
//
//		cout<<"mHChol[0] "<<mHChol[0]<<endl;
//		cout<<"mHChol[1] "<<mHChol[1]<<endl;
//		cout<<"mHChol[2] "<<mHChol[2]<<endl;
//		cout<<"mHChol[3] "<<mHChol[3]<<endl;

		MatrixMulti(mHChol, mEpsDraw,iDimOfObs, iDimOfObs, iDimOfObs, 1, mEps);
		MatrixMulti(mZ, mNextAPlus,iDimOfObs, iDimOfStates, iDimOfStates, 1, mZa);
		MatrixAdd(mZa, mEps, iDimOfObs, 1, mNextYPlus);

		/* Delete Random later just for check */

//		for(int k=0; k<iDimOfStates;k++)
//		{
//			mRandom[(2+iDimOfStates)*i+k]=mEta[k];
//		}
//
//
//		if(vN[i]==0)
//		{
//			mRandom[(2+iDimOfStates)*(i+1)+iDimOfStates+0]=mEps[0];
//			mRandom[(2+iDimOfStates)*(i+1)+iDimOfStates+1]=0;
//		}
//		else
//		{
//			for(int k=0; k<2;k++)
//			{
//				mRandom[(2+iDimOfStates)*(i+1)+iDimOfStates+k]=mEps[k];
//			}
//
//		}


		if(vN[i]==0)
		{
			iDimOfObs=2;
		}

		/* Update pointers */
		if(i!=0 )
		{
			mAPlus=mAPlus+iDimOfStates;
			mYPlus=mYPlus+iDimOfObs;
		}

		if(i<iNumOfObs-1)
		{
			mH=mH+iDimOfObs*iDimOfObs;

		}
		for(int k=0; k<iDimOfStates; k++)
		{
			mAPlus[k]=mNextAPlus[k];

		}
		for(int k=0; k<iDimOfObs; k++)
		{
			mYPlus[k]=mNextYPlus[k];
		}
	}

	//cout << "mP1[0] " <<mP1[0] <<endl;

	delete [] mQChol ;
	delete [] mHChol  ;
	delete [] mP1Chol ;
	delete [] mEpsDraw;
	delete [] mEtaDraw;
	delete [] mEps;
	delete [] mEta;
	delete [] mZa;
	delete [] mTa;
	delete [] mNextAPlus;
	delete [] mNextYPlus;

	/* Save random values */
//	string sRFile="mRandom.csv";
//	WriteOutDoubleArray(mRandom, iNumOfObs+1, (iDimOfStates+iDimOfObs), sRFile);
//	cout<<"Done"<<endl;
//	delete [] mRandom;

}

void SimulationSmootherWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam,double *mY,double *mH,unsigned int iNumOfObs, unsigned int iDimOfObs,unsigned int iDimOfStates, const gsl_rng * gsl_random_num, double* mDraw)
{
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatricesWithSeasonal(iNumOfKnots,sParam,mZ, mT, mQ, mA1, mP1);

	double * mAHat =new double[iDimOfStates*iNumOfObs];
	double * mVHat =new double[iDimOfStates*iDimOfStates*iNumOfObs];
	double * mAHatPlus =new double[iDimOfStates*iNumOfObs];
	double * mVHatPlus =new double[iDimOfStates*iDimOfStates*iNumOfObs];
	double * mAPlus =new double[iDimOfStates*iNumOfObs];
	double * mYPlus =new double[iDimOfObs*iNumOfObs];

	double * mA =new double[iNumOfObs*iDimOfStates];
	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
	double * mV =new double[iNumOfObs*iDimOfObs];
	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];

	/* APlus */

//	string sYFile="mY.csv";
//	WriteOutDoubleArray(mY, iNumOfObs, 2, sYFile);
//	string sHFile="mH.csv";
//	WriteOutDoubleArray(mH, iNumOfObs, 4, sHFile);
//	string sNFile="mN.csv";
//	WriteOutIntArray(sParam.vN, iNumOfObs, 1, sNFile);

	SSFRecursionWithSeasonal(iNumOfKnots,    mWTilde ,sParam,sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1,  iNumOfObs,  iDimOfObs,iDimOfStates,  gsl_random_num, mYPlus, mAPlus);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAPlus at " << i<< " is "<< mAPlus[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfObs; i++)
//	{
//		cout <<"mYPlus at " << i<< " is "<< mYPlus[i]<< endl;
//	}

//	for(int i=0;i<int(iNumOfObs); i++)
//	{
//		for(int j=0; j<int(iDimOfStates); j++)
//		{
//			if(isthisnan(mAPlus[i*iDimOfStates+j])==1)
//			{
//				cout<<"Obs "<< i<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2]<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2+1]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4+3]<<endl;
//			}
//		}
//
//		for(int j=0; j<int(iDimOfObs); j++)
//		{
//			if(isthisnan(mYPlus[i*iDimOfStates+j])==1)
//			{
//				cout<<"Obs "<< i<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2]<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2+1]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4+3]<<endl;
//			}
//		}
//	}

	/* AHat */
	KalmanFilterWithSeasonal(iNumOfKnots,   mWTilde ,mY, sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
	KalmanSmootherWithSeasonal( iNumOfKnots,    mWTilde ,sParam.vN,mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHat,  mVHat );


//	for(int i=0;i<int(iNumOfObs); i++)
//	{
//		for(int j=0; j<int(iDimOfStates); j++)
//		{
//			if((isthisnan(mAHat[i*iDimOfStates+j])==1)|(isthisnan(mA[i*iDimOfStates+j])==1))
//			{
//				cout<<"Obs "<< i<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2]<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2+1]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4+3]<<endl;
//			}
//		}
//		for(int j=0; j<int(iDimOfStates*iDimOfStates); j++)
//		{
//			if((isthisnan(mVHat[i*iDimOfStates*iDimOfStates+j])==1)|(isthisnan(mP[i*iDimOfStates*iDimOfStates+j])==1))
//			{
//				cout<<"Obs "<< i<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2]<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2+1]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4+3]<<endl;
//			}
//
//		}
//	}


//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mY at " << i<< " is "<< mY[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mH at " << i<< " is "<< mH[i]<< endl;
//	}
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAHat at " << i<< " is "<< mAHat[i]<< endl;
//	}

	/* AHatPlus */

	KalmanFilterWithSeasonal(iNumOfKnots,   mWTilde ,mYPlus, sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mV at " << i<< " is "<< mV[i]<< endl;
//	}

	KalmanSmootherWithSeasonal( iNumOfKnots,    mWTilde ,sParam.vN, mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHatPlus,  mVHatPlus);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAHatPlus at " << i<< " is "<< mAHatPlus[i]<< endl;
//	}

//
//	for(int i=0;i<int(iNumOfObs); i++)
//	{
//		for(int j=0; j<int(iDimOfStates); j++)
//		{
//			if((isthisnan(mAHatPlus[i*iDimOfStates+j])==1)|(isthisnan(mA[i*iDimOfStates+j])==1))
//			{
//				cout<<"Obs Plus"<< i<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2]<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2+1]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4+3]<<endl;
//			}
//		}
//		for(int j=0; j<int(iDimOfStates*iDimOfStates); j++)
//		{
//			if((isthisnan(mVHatPlus[i*iDimOfStates*iDimOfStates+j])==1)|(isthisnan(mP[i*iDimOfStates*iDimOfStates+j])==1))
//			{
//				cout<<"Obs "<< i<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2]<<endl;
//				cout<<"mAuxY  "<<sParam.mAuxY[i*2+1]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4]<<endl;
//				cout<<"mAuxH  "<<sParam.mAuxH[i*4+3]<<endl;
//			}
//		}
//	}

	double * mAHatAHatPlus =new double[iDimOfStates*iNumOfObs];
	MatrixSub(mAHat, mAHatPlus, iNumOfObs, iDimOfStates, mAHatAHatPlus);
	MatrixAdd(mAHatAHatPlus, mAPlus, iNumOfObs, iDimOfStates, mDraw);

//	cout << "mDraw[0] " << mDraw[0] << endl;
//	cout << "mDraw[iNumOfObs-1] " << mDraw[2*iNumOfObs-1] << endl;

	delete [] mZ;
	delete [] mT;
	delete [] mQ;
	delete [] mA1;
	delete [] mP1;

	delete [] mAHat;
	delete [] mVHat;
	delete [] mAHatPlus;
	delete [] mVHatPlus;
	delete [] mAPlus;
	delete [] mYPlus;

	delete [] mA;
	delete [] mP;
	delete [] mFInv;
	delete [] mV;
	delete [] mL;
	delete [] mAHatAHatPlus;
}

void DrawXandMuWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam, int iNumOfObs, const gsl_rng * gsl_random_num )
// 		DrawXandMuWithSeasonal(iNumOfKnots, mWTilde, sParam, iNumOfObs, gsl_random_num );

{
	int iDimOfObs=2;
	int iDimOfStates= iNumOfKnots-1+2;
	double * mDraw=new double[ iDimOfStates*iNumOfObs];

//	for(int i=0; i<iNumOfObs; i++)
//	{
//		sParam.mAuxY[2*i]=sParam.mAuxY[2*i]-sParam.vS[i];
//		sParam.mAuxY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS[i];
//	}
	SimulationSmootherWithSeasonal(iNumOfKnots,   mWTilde ,sParam,sParam.mAuxY,sParam.mAuxH, iNumOfObs,iDimOfObs,iDimOfStates, gsl_random_num, mDraw);

//	double * mY=new double[2*iNumOfObs];
//	double * mH=new double[4*iNumOfObs];
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		mY[2*i]=sParam.mAuxY[2*i]-sParam.vS[i];
//		mY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS[i];
//	}
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		mH[4*i]=sParam.mAuxH[4*i];
//		mH[4*i+1]=sParam.mAuxH[4*i+1];
//		mH[4*i+2]=sParam.mAuxH[4*i+2];
//		mH[4*i+3]=sParam.mAuxH[4*i+3];
//	}
//
							// for(int j=0; j<iNumOfKnots-1; j++)
							// {
								// cout <<"beta " <<mDraw[j]<<endl;
							// }

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(mWTilde,i,j)*mDraw[j];
		}
	}

	for(int j=0; j<iNumOfKnots-1;j++)
	{
		sParam.vBeta[j]=mDraw[j];
	}

	sParam.dMu[0]=mDraw[iDimOfStates-2];

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vX[i]=mDraw[iDimOfStates*i+iDimOfStates-1];
	}

	delete [] mDraw;
}

double PhiSigmaLogPosteriorWithSeasonal(const gsl_vector *v, void *ParamPointer)
{
	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double dLogitPhi= gsl_vector_get(v, 0);
	double dPhi=LogitTransformBack(dLogitPhi);
	double dLogSigma2= gsl_vector_get(v, 1);
	double dSigma2=exp(dLogSigma2);
	double dLLValue;

	double dPhiPrior=log(gsl_ran_beta_pdf ((dPhi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	double dSigmaPrior=log(gsl_ran_gamma_pdf (1/dSigma2, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));

//	cout <<"(Param->iNumOfObs)[0] in function GSL " <<  (Param->iNumOfObs)[0] << endl;
	/*LogLikelihood */
	double *dLL= new double;
	double dPhiTrue=(Param->dPhi)[0];
	double dSigma2True=(Param->dSigma2)[0];
	(Param->dPhi)[0]=dPhi;
	(Param->dSigma2)[0]=dSigma2;

	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  dLL );

	(Param->dSigma2)[0]=dSigma2True;
	(Param->dPhi)[0]=dPhiTrue;
	dLLValue=dLL[0];
	delete dLL;
//	cout <<"Value dPhi " <<dPhi << endl;
//	cout <<"Value dSigma " << dSigma << endl;
//	cout <<"Value " <<-(dLLValue+dPhiPrior+dSigmaPrior)/((Param->iNumOfObs)[0]) <<endl;
	return  -(dLLValue+dPhiPrior+dSigmaPrior)/((Param->iNumOfObs)[0]);
}

void PhiSigmaLogPosteriorDerivativeWithSeasonal(const gsl_vector *v, void *ParamPointer, gsl_vector *df)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= gsl_vector_get(v, 0);
	double y= gsl_vector_get(v, 1);

	double hforx=pow(dEps,1.0/3)*fmax( fabs(x),1);
	double hfory=pow(dEps,1.0/3)*fmax( fabs(y),1);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);


//	cout<<"DERIV phi0 "<< phi0 << endl;
//	cout<<"DERIV phi "<< phi << endl;
//	cout<<"DERIV phi1 "<< phi1 << endl;
//	cout<<"DERIV sigma0 "<< sigma0 << endl;
//	cout<<"DERIV sigma "<< sigma << endl;
//	cout<<"DERIV sigma1 "<< sigma1 << endl;


	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

	/* derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);


//	cout<<"DERIV fxh0 "<< fxh0 << endl;
//	cout<<"DERIV fxh1 "<< fxh1 << endl;
//	cout<<"DERIV hhforx "<< hhforx << endl;
//	cout<<"DERIV fxh1-fxh0 "<< fxh1-fxh0 << endl;
//	cout<<"DERIV (fxh1-fxh0)/hhforx "<< (fxh1-fxh0)/(2*hhforx) << endl;


	gsl_vector_set(df, 0, (fxh1-fxh0)/(2*hhforx));


	/* derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);


//	cout<<"DERIV fyh0 "<< fyh0 << endl;
//	cout<<"DERIV fyh1 "<< fyh1 << endl;
//	cout<<"DERIV hhfory "<< hhfory << endl;
//	cout<<"DERIV fyh1-fyh0 "<< fyh1-fyh0 << endl;
//	cout<<"DERIV (fyh1-fyh0)/hhfory "<< (fyh1-fyh0)/(2*hhfory) << endl;

	gsl_vector_set(df, 1,(fyh1-fyh0)/(2*hhfory));

	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;
//	cout<<"XXXXXXXXXXXXXXXXXXXXXX DERIV END XXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
}

void PhiSigmaLogPosteriorANDDerivativeWithSeasonal(const gsl_vector *v, void *ParamPointer,
        double *f, gsl_vector *df)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= gsl_vector_get(v, 0);
	double y= gsl_vector_get(v, 1);

	double hforx=pow(dEps,1.0/3)*fmax( fabs(x),1);
	double hfory=pow(dEps,1.0/3)*fmax( fabs(y),1);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);

	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

	/* function evaluation */
	double fx=0;
	fx=fx+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fx=fx+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fx=-(fx+LL[0])/((Param->iNumOfObs)[0]);

	f[0]=fx;

	/* derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);

	gsl_vector_set(df, 0, (fxh1-fxh0)/(2*hhforx));


	/* derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);
	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);

	gsl_vector_set(df, 1,(fyh1-fyh0)/(2*hhfory));

	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;

}

void CalculatingHessianWithSeasonal(const gsl_vector *v, void *ParamPointer,double * mHessian)
{
	double dEps=2.22044604925031308e-16;

	struct AllParam * Param = (struct AllParam *)ParamPointer;
	double x= gsl_vector_get(v, 0);
	double y= gsl_vector_get(v, 1);


	double hforx=pow(dEps,1.0/4)*fmax( (double) fabs((double) x),1.0);
	double hfory=pow(dEps,1.0/4)*fmax( (double) fabs( (double) y),1.0);

	double xh1=x+hforx;
	double xh0=x-hforx;

	double yh1=y+hfory;
	double yh0=y-hfory;

	double hhforx=xh1-xh0;
	double hhfory=yh1-yh0;


	double phi1=LogitTransformBack(xh1);
	double phi0=LogitTransformBack(xh0);
	double phi=LogitTransformBack(x);
	double sigma1=exp(yh1);
	double sigma0=exp(yh0);
	double sigma=exp(y);

	double phiTrue=(Param->dPhi)[0];
	double sigmaTrue=(Param->dSigma2)[0];

	double* LL= new double;

//	cout<<"HESSIAN phi0 "<< phi0 << endl;
//	cout<<"HESSIAN phi "<< phi << endl;
//	cout<<"HESSIAN phi1 "<< phi1 << endl;
//	cout<<"HESSIAN sigma0 "<< sigma0 << endl;
//	cout<<"HESSIAN sigma "<< sigma << endl;
//	cout<<"HESSIAN sigma1 "<< sigma1 << endl;

	/* function evaluation */
	double fx=0;
	fx=fx+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fx=fx+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fx=-(fx+LL[0])/((Param->iNumOfObs)[0]);

	/* 2nd derivative with respect to x */
	/* forward */
	double fxh1=0;
	fxh1=fxh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1=fxh1+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh1=-(fxh1+LL[0])/((Param->iNumOfObs)[0]);

	/* backward */
	double fxh0=0;
	fxh0=fxh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0=fxh0+log(gsl_ran_gamma_pdf (1/sigma, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh0=-(fxh0+LL[0])/((Param->iNumOfObs)[0]);

//	cout<<"HESSIAN fxh0 "<< fxh0 << endl;
//	cout<<"HESSIAN fxh1 "<< fxh1 << endl;
//	cout<<"HESSIAN fx "<< fx << endl;
//	cout<<"HESSIAN hhforx "<< hhforx << endl;
//	cout<<"HESSIAN fxh1-2*fx+fxh0 "<< fxh1-2*fx+fxh0<< endl;
//	cout<<"HESSIAN (fxh1-2*fx+fxh0)/pow(hhforx,2)"<< (fxh1-2*fx+fxh0)/pow(hhforx,2)<< endl;

	mHessian[0]= (Param->iNumOfObs)[0]*((fxh1-2*fx+fxh0)/pow(hhforx,2));

	/* 2nd derivative with respect to y */
	/* forward */
	double fyh1=0;
	fyh1=fyh1+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh1=fyh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma1;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fyh1=-(fyh1+LL[0])/((Param->iNumOfObs)[0]);

	/* backward */
	double fyh0=0;
	fyh0=fyh0+log(gsl_ran_beta_pdf ((phi+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fyh0=fyh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi;
	(Param->dSigma2)[0]=sigma0;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fyh0=-(fyh0+LL[0])/((Param->iNumOfObs)[0]);

	mHessian[3]= (Param->iNumOfObs)[0]*((fyh1-2*fx+fyh0)/pow(hhfory,2));


	/* corss derivative with respect to x and y */
	/* x forward  y forward*/
	double fxh1yh1=0;
	fxh1yh1=fxh1yh1+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1yh1=fxh1yh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma1;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh1yh1=-(fxh1yh1+LL[0])/((Param->iNumOfObs)[0]);
	/* x forward  y backward*/
	double fxh1yh0=0;
	fxh1yh0=fxh1yh0+log(gsl_ran_beta_pdf ((phi1+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh1yh0=fxh1yh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi1;
	(Param->dSigma2)[0]=sigma0;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh1yh0=-(fxh1yh0+LL[0])/((Param->iNumOfObs)[0]);
	/* x backward  y forward*/
	double fxh0yh1=0;
	fxh0yh1=fxh0yh1+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0yh1=fxh0yh1+log(gsl_ran_gamma_pdf (1/sigma1, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma1;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh0yh1=-(fxh0yh1+LL[0])/((Param->iNumOfObs)[0]);
	/* x backward  y backward*/
	double fxh0yh0=0;
	fxh0yh0=fxh0yh0+log(gsl_ran_beta_pdf ((phi0+1)/2, (Param->dPriorPhiA)[0], (Param->dPriorPhiB)[0]));
	fxh0yh0=fxh0yh0+log(gsl_ran_gamma_pdf (1/sigma0, (Param->dPriorSigmaA)[0], (Param->dPriorSigmaB)[0]));
	(Param->dPhi)[0]=phi0;
	(Param->dSigma2)[0]=sigma0;
	CalculateLLWithSeasonal((Param->iNumOfKnots)[0], (Param->W), Param[0],(Param->iNumOfObs)[0], 2, (Param->iNumOfKnots)[0]-1+2,  LL );
	fxh0yh0=-(fxh0yh0+LL[0])/((Param->iNumOfObs)[0]);

	mHessian[1]=(Param->iNumOfObs)[0]*((fxh1yh1-fxh1yh0-fxh0yh1+fxh0yh0)/(4*hhforx*hhfory));
	mHessian[2]=(Param->iNumOfObs)[0]*((fxh1yh1-fxh1yh0-fxh0yh1+fxh0yh0)/(4*hhforx*hhfory));


	(Param->dPhi)[0]=phiTrue;
	(Param->dSigma2)[0]=sigmaTrue;
	delete LL;


}

void DrawPhiSigmaLaplaceApproxWithSeasonal( struct AllParam sParam ,  int iNumOfObs, const gsl_rng * gsl_random_num )
{


	/* Setting up the optimizer */
	int iter = 0;
	int max_iter=100;
	int status;

//	cout << "Num Of obs : " << iNumOfObs << endl;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf F;


	/* setting up the objective function */
	F.n = 2;
	F.f =&PhiSigmaLogPosteriorWithSeasonal;
	F.df = &PhiSigmaLogPosteriorDerivativeWithSeasonal;
	F.fdf = &PhiSigmaLogPosteriorANDDerivativeWithSeasonal;
	F.params = &sParam;

	/* Starting point*/
	x = gsl_vector_alloc (2);
	gsl_vector_set (x,0, LogitTransform(sParam.dPhi[0]));
	gsl_vector_set (x,1,log(sParam.dSigma2[0]));

	/* Setting up the minimizer */
	T =  gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, 2);

	gsl_multimin_fdfminimizer_set (s, &F, x, 0.000001, 0.000001);

	/* Iterating the optimizer */
	do
	{
	      iter++;
	      status = gsl_multimin_fdfminimizer_iterate (s);

//	      cout<<iter<<endl;
	      if (status)
	        break;

	      status = gsl_multimin_test_gradient (s->gradient, 1e-2);

//	      if (status == GSL_SUCCESS)
//	      printf ("Minimum found at:\n");
//
//	      printf ("%5d %.5f %10.5f\n", iter,
//	              gsl_vector_get (s->x, 0),
//	              s->f);

	  }
	  while (status == GSL_CONTINUE && iter < max_iter);


//	  /* Mean */
	  double* mMeanLogit=new double[2];
	  mMeanLogit[0]=gsl_vector_get (s->x, 0);
	  mMeanLogit[1]=gsl_vector_get(s->x,1);


	  double* mInfoMatrix=new double[4];
//	  double* mInfoMatrixJacobian=new double[4];
	  gsl_vector *atx;
	  atx=gsl_vector_alloc (2);
	  gsl_vector_set (atx,0, gsl_vector_get (s->x, 0));
	  gsl_vector_set (atx,1, gsl_vector_get (s->x, 1));
	  CalculatingHessianWithSeasonal(atx, &sParam, mInfoMatrix);


	  double dDet=mInfoMatrix[0]*mInfoMatrix[3]-mInfoMatrix[1]*mInfoMatrix[2];
	  if((dDet<0) | (isthisnan(dDet)==1))
	  {
		  mInfoMatrix[1]=0;
		  mInfoMatrix[2]=0;
	  }

	  if((mInfoMatrix[0]<0) |( isthisnan(mInfoMatrix[0])==1) )
	  {
		  mInfoMatrix[0]=0.01;
	  }
	  if( (mInfoMatrix[3]<0) | (isthisnan(mInfoMatrix[3])==1) )
	  {
		  mInfoMatrix[3]=0.01;
	  }
	  /* Inverse info matrix */
	  double* mInvInfoMatrix=new double[4];
	  MatirxInv2x2(mInfoMatrix , 2,mInvInfoMatrix);
//	  if(dDet<0)
//	  {
//	 		  mInvInfoMatrix[1]=0.075;
//	 		  mInvInfoMatrix[2]=0.075;
//	  }
//	  for(int k=0; k<4; k++)
//	  {
//	  		  cout<<"mInvInfoMatrix at "<< k <<" is "<<mInvInfoMatrix[k]<<endl;
//	  }

	  /* Cholesky Inverse info matrix  */
	  double* mCholInvInfoMatrix=new double[4];
//	  cout<<"Phi Sigma Laplace chol start"<<endl;
	  Cholesky2x2(mInvInfoMatrix, 2, mCholInvInfoMatrix);
//	  cout<<"Phi Sigma Laplace chol end"<<endl;
//	  for(int k=0; k<4; k++)
//	  {
//	 	  		  cout<<"mCholInvInfoMatrix at "<< k <<" is "<<mCholInvInfoMatrix[k]<<endl;
//	  }


	  /* Draw new Phi and Sigma from t density */
	  double* mRandom=new double[2];
	  double* mNewLogit=new double[2];
	  unsigned int  iDf=5;
	  double dSqrtW=sqrt((double) iDf/gsl_ran_chisq (gsl_random_num, iDf));

	  for (int k=0;k<2;k++)
	  {
		  mRandom[k]=dSqrtW*gsl_ran_gaussian(gsl_random_num,1);

	  }
	  MatrixMulti(mCholInvInfoMatrix, mRandom,2, 2, 2, 1, mNewLogit);
	  MatrixAdd(mMeanLogit, mNewLogit, 2, 1, mNewLogit);



	  gsl_vector * newx;
	  newx=gsl_vector_alloc(2);
	  gsl_vector_set (newx,0, mNewLogit[0]);
	  gsl_vector_set (newx,1, mNewLogit[1]);

	  gsl_vector * oldx;
	  oldx=gsl_vector_alloc(2);
	  gsl_vector_set (oldx,0, LogitTransform(sParam.dPhi[0]));
	  gsl_vector_set (oldx,1, log(sParam.dSigma2[0]));
	  double * mOldLogit=new double[2];
	  mOldLogit[0]=LogitTransform(sParam.dPhi[0]);
	  mOldLogit[1]=log(sParam.dSigma2[0]);


	  /*Calculate Acceptance rate */
	  double dLogU=log(gsl_rng_uniform(gsl_random_num));

//	  double sigma_dPhi=sqrt(mInvInfoMatrix[0]);
//	  double sigma_dSigma=sqrt(mInvInfoMatrix[3]);
//	  double rho_dPhiSigma=mInvInfoMatrix[1]/(sigma_dPhi*sigma_dSigma);


//	  cout << "new Phi " <<LogitTransformBack( mNewLogit[0]) <<endl;
//	  cout << "new Sigma " <<exp( mNewLogit[1]) <<endl;
//	  cout << "LL at new "<<-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs <<endl;
//	  cout << "LL at old "<< +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs <<endl;
//	  cout << "proposal at new "<< log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)- mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma))<<endl;
//	  cout << "proposal at old "<< log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma)) <<endl;
//	  cout << "acceptance "<< -PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (oldx, 0)- mMeanLogit[0], gsl_vector_get (oldx, 1)- mMeanLogit[1], sigma_dPhi,sigma_dSigma, rho_dPhiSigma))
//							 		 +PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(gsl_ran_bivariate_gaussian_pdf (gsl_vector_get (newx, 0)- mMeanLogit[0], gsl_vector_get (newx, 1)- mMeanLogit[1], sigma_dPhi, sigma_dSigma, rho_dPhiSigma)) <<endl;



	  double dNewLog=-PhiSigmaLogPosteriorWithSeasonal(newx, &sParam )*iNumOfObs+ log(BivariateStudentTDensity(mOldLogit,  mMeanLogit, mInvInfoMatrix , iDf));
	  double dOldLog= -PhiSigmaLogPosteriorWithSeasonal(oldx, &sParam)*iNumOfObs+ log(BivariateStudentTDensity(mNewLogit,mMeanLogit, mInvInfoMatrix , iDf));
	  double dAccept=fmin(dNewLog-dOldLog,0);

//	  cout  <<"dLogU " << dLogU<<" dAccept " << dAccept <<endl;


//	  cout <<"new phi "<<LogitTransformBack(mNewLogit[0])<<endl;
//	  cout <<"old phi "<<LogitTransformBack(mOldLogit[0])<<endl;
//	  cout <<"new phi log"<<-PhiSigmaLogPosteriorGSL(newx, &sParam )*iNumOfObs+ log(BivariateStudentTDensity(mOldLogit,  mMeanLogit, mInvInfoMatrix , iDf))<<endl;
//	  cout <<"old phi log "<<+PhiSigmaLogPosteriorGSL(oldx, &sParam)*iNumOfObs- log(BivariateStudentTDensity(mNewLogit,mMeanLogit, mInvInfoMatrix , iDf))<<endl;
//	  cout <<"Accept "<<dAccept<<endl;
	  if((dLogU<=dAccept) & (LogitTransformBack(mNewLogit[0])!=1) & (isthisnan(dNewLog)!=1)  )
	  {

//		  printf ("new phi: %1.12f  \n", LogitTransformBack(mNewLogit[0]) );
		  sParam.dPhi[0]=LogitTransformBack(mNewLogit[0]);
		  sParam.dSigma2[0]=exp(mNewLogit[1]);

	  }
//	  cout  <<"Accepted Phi " << sParam.dPhi[0] <<" Sigma "<< sParam.dSigma2[0] << endl;






	  //cout << sParam.dGamma[0] << endl;
	  gsl_multimin_fdfminimizer_free (s);
	  gsl_vector_free (x);

	  delete [] mInfoMatrix;
	  delete [] mInvInfoMatrix;
	  delete [] mCholInvInfoMatrix;
	  delete [] mRandom;
	  delete [] mNewLogit;
	  delete [] mMeanLogit;



}

void PhiSigmaLogPosteriorWithSeasonal(int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam, int iNumOfObs, int iDimOfObs,
		 int iDimOfStates, double * dLogPost)
{
	/*Log Priors */
	/* Phi Prior */
	dLogPost[0]=log(gsl_ran_beta_pdf ((sParam.dPhi[0]+1)/2, sParam.dPriorPhiA[0], sParam.dPriorPhiB[0]));

//	cout<<"Phi prior "<<dLogPost[0]<<endl;
	/* Sigma2 Prior */
	dLogPost[0]=dLogPost[0]+log(gsl_ran_gamma_pdf (1/sParam.dSigma2[0], sParam.dPriorSigmaA[0], sParam.dPriorSigmaB[0]));
//	cout<<"Sigma prior "<<dLogPost[0]<<endl;
	/*LogLikelihood */
	double *dLL= new double;
	CalculateLLWithSeasonal( iNumOfKnots,   mWTilde , sParam,iNumOfObs, iDimOfObs, iDimOfStates,  dLL );

	dLogPost[0]=dLogPost[0]+dLL[0];
//	cout<<"dLL "<<dLL[0]<<endl;
//	cout<<"full posterior "<<dLogPost[0]<<endl;
	delete dLL;

}

void DrawPhiSigmaAdaptiveRWWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam,int iNumOfObs,  int iNumOfIter,
		const gsl_rng * gsl_random_num, double* mCovar,  double* mSum )
{
	/* Proposal */
	double dOmega1=0.05;
	double dU=gsl_rng_uniform(gsl_random_num);
//	double dNewPhiLogit;
	double* mChol=new double[4];
	double* mSigma=new double[4];
	double* mSumSum=new double[4];
	double* mParamParam=new double[4];
	double *mSumTrans= new double[2];
	double *mLogParam=new double[2];
	double *mLogParamTrans=new double[2];
	double *mNewLogParam=new double[2];
	double *mRandom=new double[2];
	int iFlag=1;



	if(iNumOfIter<=1)
	{
//		cout<<"first two iter"<<endl;
		mChol[0]=0.1/sqrt(2);
		mChol[1]=0;
		mChol[2]=0;
		mChol[3]=0.1/sqrt(2);
	}
	else if(iNumOfIter==2)
	{
//		cout<<"third iter"<<endl;
		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
		MatrixScalarMulti(mSumSum,0.5*0.5,2, 2, mSumSum);
		MatrixScalarMulti(mCovar,0.5,2, 2, mCovar);
		MatrixSub(mCovar,mSumSum, 2,2, mSigma);
		iFlag=Cholesky2x2(mSigma,2, mChol);
//		cout<<"mSigma[0] "<<mSigma[0]<<endl;
//		cout<<"mSigma[1] "<<mSigma[1]<<endl;
//		cout<<"mSigma[2] "<<mSigma[2]<<endl;
//		cout<<"mSigma[3] "<<mSigma[3]<<endl;
	}
	else
	{
//		cout<<"later iter"<<endl;
		iFlag=Cholesky2x2(mCovar,2, mChol);

	}

	mLogParam[0]=LogitTransform(sParam.dPhi[0]);
	mLogParam[1]=log(sParam.dSigma2[0]);

	if(dU<=dOmega1)
	{
		/* Static proposal */
		mChol[0]=0.1/sqrt(2);
		mChol[1]=0;
		mChol[2]=0;
		mChol[3]=0.1/sqrt(2);

		mRandom[0]=gsl_ran_gaussian(gsl_random_num,1);
		mRandom[1]=gsl_ran_gaussian(gsl_random_num,1);

//		cout<<"mRandom[0] "<<mRandom[0]<<endl;
//		cout<<"mRandom[1] "<<mRandom[1]<<endl;

		MatrixMulti(mChol, mRandom,2, 2, 2, 1, mNewLogParam);
		MatrixAdd(mNewLogParam, mLogParam, 2, 1, mNewLogParam);
	}
	else
	{
		mRandom[0]=gsl_ran_gaussian(gsl_random_num,1);
		mRandom[1]=gsl_ran_gaussian(gsl_random_num,1);
//		cout<<"mRandom[0] "<<mRandom[0]<<endl;
//		cout<<"mRandom[1] "<<mRandom[1]<<endl;
//		cout<<"mLogParam[0] "<<mLogParam[0]<<endl;
//		cout<<"mLogParam[1] "<<mLogParam[1]<<endl;

		if(iFlag==0 )
		{
			mChol[0]=0.1/sqrt(2);
			mChol[1]=0;
			mChol[2]=0;
			mChol[3]=0.1/sqrt(2);
		}
		else
		{
			MatrixScalarMulti(mChol,2.38/sqrt(2),2, 2, mChol);
		}
//		cout<<"iFlag "<<iFlag<<endl;
//		cout<<"mChol[0] "<<mChol[0]<<endl;
//		cout<<"mChol[1] "<<mChol[1]<<endl;
//		cout<<"mChol[2] "<<mChol[2]<<endl;
//		cout<<"mChol[3] "<<mChol[3]<<endl;


//		mChol[0] =0.001;
//		mChol[1] =0;
//		mChol[2] =0;
//		mChol[3] =0.0005;

		MatrixMulti(mChol, mRandom,2, 2, 2, 1, mNewLogParam);
		MatrixAdd(mNewLogParam, mLogParam, 2, 1, mNewLogParam);
	}



	/* Acceptance Rate */

	double dLogU=log(gsl_rng_uniform(gsl_random_num));

	double* dNewLogPosterior= new double;
	double* dOldLogPosterior= new double;


	sParam.dPhi[0]=LogitTransformBack(mNewLogParam[0]);
	sParam.dSigma2[0]=exp(mNewLogParam[1]);

//	cout<<"XXXXX new XXXXXXXX"<<endl;
//	cout<<" new Phi " << sParam.dPhi[0]<< endl;
//	cout<<" new Sigma " << sParam.dSigma[0]<< endl;
	PhiSigmaLogPosteriorWithSeasonal(iNumOfKnots,   mWTilde,sParam, iNumOfObs, 2, iNumOfKnots-1+2,dNewLogPosterior);

	sParam.dPhi[0]=LogitTransformBack(mLogParam[0]);
	sParam.dSigma2[0]=exp(mLogParam[1]);

//	cout<<"XXXXX old XXXXXXXX"<<endl;
//	cout<<" old Phi " << sParam.dPhi[0]<< endl;
//	cout<<" old Sigma2 " << sParam.dSigma2[0]<< endl;
	PhiSigmaLogPosteriorWithSeasonal(iNumOfKnots,  mWTilde,sParam, iNumOfObs, 2, iNumOfKnots-1+2,dOldLogPosterior);

    double dAccept=fmin(dNewLogPosterior[0]-dOldLogPosterior[0],0);


//    cout<< "old Phi " << LogitTransformBack(mLogParam[0]) << " new Phi " << LogitTransformBack(mNewLogParam[0]) << endl;
//    cout<< "old Sigma2 " << exp(mLogParam[1]) << " new Sigma2 " << exp(mNewLogParam[1])<< endl;
//    cout<< "old " << dOldLogPosterior[0]<< " new " << dNewLogPosterior[0]<< endl;
//    cout<< "dLogU " << dLogU<<" dAccept " << dAccept<<  endl;


	if((dLogU<=dAccept) & (LogitTransformBack(mNewLogParam[0])!=1) & (isthisnan(dNewLogPosterior[0])!=1) )
	{
		sParam.dPhi[0]=LogitTransformBack(mNewLogParam[0]);
		sParam.dSigma2[0]=exp(mNewLogParam[1]);
//		 cout<< "XXXXXXXX Accepted XXXXXXXXX"<< endl;
	}
//	cout<< "XXXXXXXXXXXXXXXXX"<< endl;

//	cout<<"final Phi " << sParam.dPhi[0]<< endl;
//	cout<<"final Sigma2 " << sParam.dSigma2[0]<< endl;

	mLogParam[0]=LogitTransform(sParam.dPhi[0]);
	mLogParam[1]=log(sParam.dSigma2[0]);
	double dNumOfIter=(double) iNumOfIter;
	/* Covariance */
	if(iNumOfIter<=1)
	{
//		cout<<"first two iter END"<<endl;
//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;
		MatrixTrans( mLogParam, 2, 1, mLogParamTrans);
		MatrixMulti(mLogParam,mLogParam,2, 1,1,2, mParamParam);
		MatrixAdd(mCovar, mParamParam, 2,2, mCovar);
		MatrixAdd(mSum, mLogParam, 2,1, mSum);
//		cout<<"mSum[0] "<<mSum[0]<<endl;
//		cout<<"mSum[1] "<<mSum[1]<<endl;
//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

	}
	else if(iNumOfIter==2)
	{
//		cout<<"third iter END"<<endl;
		MatrixScalarMulti(mSigma,(dNumOfIter-1)/dNumOfIter,2,2,mCovar);

//		cout<<"mSigma[0] "<<mSigma[0]<<endl;
//		cout<<"mSigma[1] "<<mSigma[1]<<endl;
//		cout<<"mSigma[2] "<<mSigma[2]<<endl;
//		cout<<"mSigma[3] "<<mSigma[3]<<endl;
//
//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
//		cout<<"SumSum[0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum[1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum[2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum[3] "<<mSumSum[3]<<endl;
		MatrixScalarMulti(mSumSum,1/(dNumOfIter* dNumOfIter),2, 2, mSumSum);
//		cout<<"SumSum divided [0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum divided [1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum divided [2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum divided [3] "<<mSumSum[3]<<endl;
		MatrixAdd(mCovar, mSumSum, 2,2, mCovar);

//		cout<<"mSum[0] "<<mSum[0]<<endl;
//		cout<<"mSum[1] "<<mSum[1]<<endl;
//		cout<<"mCovar+SumSum[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum[3] "<<mCovar[3]<<endl;

		MatrixTrans( mLogParam, 2, 1, mLogParamTrans);
		MatrixMulti(mLogParam,mLogParamTrans,2, 1,1,2, mParamParam);
		MatrixScalarMulti(mParamParam,1/dNumOfIter,2, 2, mParamParam);
		MatrixAdd(mCovar, mParamParam, 2,2, mCovar);

//		cout<<"mParamParam[0] "<<mParamParam[0]<<endl;
//		cout<<"mParamParam[1] "<<mParamParam[1]<<endl;
//		cout<<"mParamParam[2] "<<mParamParam[2]<<endl;
//		cout<<"mParamParam[3] "<<mParamParam[3]<<endl;
//
//		cout<<"mCovar+SumSum+ParamParam[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[3] "<<mCovar[3]<<endl;
//
//
//		cout<<"mSum[0] "<<mSum[0]<<endl;
//		cout<<"mSum[1] "<<mSum[1]<<endl;
//		cout<<"mLogParam[0] "<<mLogParam[0]<<endl;
//		cout<<"mLogParam[1] "<<mLogParam[1]<<endl;
		mSum[0]=mSum[0]+ mLogParam[0];
		mSum[1]=mSum[1]+ mLogParam[1];

//		cout<<"mSum[0]+log"<<mSum[0]<<endl;
//		cout<<"mSum[1]+log"<<mSum[1]<<endl;

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
//		cout<<"SumSum[0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum[1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum[2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum[3] "<<mSumSum[3]<<endl;
		MatrixScalarMulti(mSumSum,1/(dNumOfIter*(dNumOfIter+1)),2, 2, mSumSum);
		MatrixSub(mCovar, mSumSum, 2,2, mCovar);
//		cout<<"SumSum divided [0] "<<mSumSum[0]<<endl;
//		cout<<"SumSum divided [1] "<<mSumSum[1]<<endl;
//		cout<<"SumSum divided [2] "<<mSumSum[2]<<endl;
//		cout<<"SumSum divided [3] "<<mSumSum[3]<<endl;

//		cout<<"mCovar+SumSum+ParamParam+last[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[3] "<<mCovar[3]<<endl;


	}
	else
	{
//		cout<<"later iter END"<<endl;

//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

		MatrixScalarMulti(mCovar,(dNumOfIter-1)/dNumOfIter,2,2,mCovar);

//		cout<<"mCovar[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar[3] "<<mCovar[3]<<endl;

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
		MatrixScalarMulti(mSumSum,1/(dNumOfIter* dNumOfIter),2, 2, mSumSum);
		MatrixAdd(mCovar, mSumSum, 2,2, mCovar);

//		cout<<"mCovar+SumSum[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum[3] "<<mCovar[3]<<endl;

		MatrixTrans( mLogParam, 2, 1, mLogParamTrans);
//		cout<<"mLogParam[0] "<<mLogParam[0]<<endl;
//		cout<<"mLogParam[1] "<<mLogParam[1]<<endl;

		MatrixMulti(mLogParam,mLogParamTrans,2, 1,1,2, mParamParam);

//		cout<<"mParamParam[0] "<<mParamParam[0]<<endl;
//		cout<<"mParamParam[1] "<<mParamParam[1]<<endl;
//		cout<<"mParamParam[2] "<<mParamParam[2]<<endl;
//		cout<<"mParamParam[3] "<<mParamParam[3]<<endl;

		MatrixScalarMulti(mParamParam,1/( dNumOfIter),2, 2, mParamParam);
		MatrixAdd(mCovar, mParamParam, 2,2, mCovar);

//		cout<<"mCovar+SumSum+ParamParam[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam[3] "<<mCovar[3]<<endl;

		mSum[0]=mSum[0]+LogitTransform(sParam.dPhi[0]);
		mSum[1]=mSum[1]+log(sParam.dSigma2[0]);

		MatrixTrans( mSum, 2, 1, mSumTrans);
		MatrixMulti(mSum,mSumTrans,2, 1,1,2, mSumSum);
		MatrixScalarMulti(mSumSum,1/( dNumOfIter* (dNumOfIter+1)),2, 2, mSumSum);
		MatrixSub(mCovar, mSumSum, 2,2, mCovar);
//		cout<<"mCovar+SumSum+ParamParam+last[0] "<<mCovar[0]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[1] "<<mCovar[1]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[2] "<<mCovar[2]<<endl;
//		cout<<"mCovar+SumSum+ParamParam+last[3] "<<mCovar[3]<<endl;


	}

	delete [] mChol;
	delete [] mSigma;
	delete [] mSumSum;
	delete [] mParamParam;
	delete [] mSumTrans;
	delete [] mLogParam;
	delete [] mLogParamTrans;
	delete [] mNewLogParam;
	delete [] mRandom;
	delete dOldLogPosterior;
	delete dNewLogPosterior;

}

//void TestSSF(unsigned int iNumOfObs)
//{
//
//	/* Initialize parameters and latent variables */
//	const gsl_rng_type * T;
//	gsl_rng * gsl_random_num;
//	T=gsl_rng_default;
//	gsl_random_num=gsl_rng_alloc(T);
//	gsl_rng_set(gsl_random_num, 100);
//
////	double* vZeroSkellamDens = new double[iNumOfObs];
////	double* vSkellamDens = new double[iNumOfObs];
////	double* vIndicator = new double[iNumOfObs];
//
//
//	/* Initialize Parameters*/
//	AllParam sParam;
//	sParam.dMu=new double;
//	sParam.dPhi=new double;
//	sParam.dSigma2=new double;
//	sParam.dGamma=new double;
//	sParam.dPriorGammaA=new double;
//	sParam.dPriorGammaB=new double;
//	sParam.dPriorMuMean=new double;
//	sParam.dPriorMuSigma2=new double;
//	sParam.dPriorPhiA=new double;
//	sParam.dPriorPhiB=new double;
//	sParam.dPriorSigmaA=new double;
//	sParam.dPriorSigmaB=new double;
//
//	sParam.vX=new double [iNumOfObs];
//	sParam.vN = new  int  [iNumOfObs];
//	sParam.vTau1= new double[iNumOfObs];
//	sParam.vTau2= new double[iNumOfObs];
//	sParam.mAuxY=new double[2*iNumOfObs];
//	sParam.mAuxH=new double[4*iNumOfObs];
//
//
//
//
//	sParam.dMu[0]=1;
//	sParam.dPhi[0]=0.99999;
//	sParam.dSigma2[0]=0.5;
//	sParam.dGamma[0]=0;
//
//	sParam.dPriorGammaA[0]=1.7;
//	sParam.dPriorGammaB[0]=10;
//	sParam.dPriorMuMean[0]=0;
//	sParam.dPriorMuSigma2[0]=1;
//	sParam.dPriorPhiA[0]=20;
//	sParam.dPriorPhiB[0]=1.5;
//	sParam.dPriorSigmaA[0]=2.5;
//	sParam.dPriorSigmaB[0]=1/0.025;
//
//
//	/* mY */
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		if ( i % 2 == 0 )
//		{
//			sParam.mAuxY[i*2]=0.3;
//			sParam.mAuxY[i*2+1]=0.5;
//		}
//		else
//		{
//			sParam.mAuxY[i*2]=0.4;
//			sParam.mAuxY[i*2+1]=0.6;
//		}
//	}
//
//	string sYFile="mY.csv";
//	WriteOutDoubleArray(sParam.mAuxY, iNumOfObs, 2, sYFile);
//
//	/* mH */
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		if ( i % 4== 0 )
//		{
//			sParam.mAuxH[i*4]=2;
//			sParam.mAuxH[i*4+1]=0;
//			sParam.mAuxH[i*4+2]=0;
//			sParam.mAuxH[i*4+3]=2;
//		}
//		else
//		{
//			sParam.mAuxH[i*4]=1;
//			sParam.mAuxH[i*4+1]=0;
//			sParam.mAuxH[i*4+2]=0;
//			sParam.mAuxH[i*4+3]=1;
//		}
//
//	}
//
//	string sHFile="mH.csv";
//	WriteOutDoubleArray(sParam.mAuxH, iNumOfObs, 4, sHFile);
//
//	/* mN */
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		if ( i % 8 == 0 )
//		{
//			sParam.vN[i]=0;
//
//		}
//		else
//		{
//			sParam.vN[i]=1;
//
//		}
//	}
//
//	string sNFile="mN.csv";
//	WriteOutIntArray(sParam.vN, iNumOfObs, 1, sNFile);
//	unsigned int iDimOfObs=2;
//	unsigned int iDimOfStates=2;
//
//	double * mZ =new double[iDimOfObs*iDimOfStates];
//	double * mT =new double[iDimOfStates*iDimOfStates];
//	double * mQ =new double[iDimOfStates*iDimOfStates];
//	double * mA1 =new double[iDimOfStates];
//	double * mP1 =new double[iDimOfStates*iDimOfStates];
//	SystemMatrices(sParam,mZ, mT, mQ, mA1, mP1);
//
//
//
//	double * mAHat =new double[iDimOfStates*iNumOfObs];
//	double * mVHat =new double[iDimOfStates*iDimOfStates*iNumOfObs];
//	double * mA =new double[iNumOfObs*iDimOfStates];
//	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
//	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
//	double * mV =new double[iNumOfObs*iDimOfObs];
//	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];
//
//
//
//
//	/* AHat */
//	KalmanFilter(sParam.mAuxY, sParam.vN, sParam.mAuxH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
//	KalmanSmoother( sParam.vN,mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHat,  mVHat );
//
//
//	string sAFile="mA.csv";
//	WriteOutDoubleArray(mA, iNumOfObs, 2, sAFile);
//	string sPFile="mP.csv";
//	WriteOutDoubleArray(mP, iNumOfObs, 4, sPFile);
//	string sAHatFile="mAHat.csv";
//	WriteOutDoubleArray(mAHat, iNumOfObs, 2, sAHatFile);
//	string sVHatFile="mVHat.csv";
//	WriteOutDoubleArray(mVHat, iNumOfObs, 4, sVHatFile);
//
//	/* LL */
//	double *dLL= new double;
//	CalculateLL( sParam,iNumOfObs, iDimOfObs, iDimOfStates,  dLL );
//	string sLLFile="mLL.csv";
//	WriteOutDoubleArray(dLL, 1, 1, sLLFile);
//	cout<<"LL "<<dLL[0]<<endl;
//
//	double * mDraw =new double[iDimOfStates*iNumOfObs];
//	SimulationSmoother(sParam,sParam.mAuxY,sParam.mAuxH, iNumOfObs, iDimOfObs,iDimOfStates,  gsl_random_num,  mDraw);
//	string sDrawFile="mDraw.csv";
//	WriteOutDoubleArray(mDraw, iNumOfObs, 2, sDrawFile);
//
//
//	delete [] mDraw;
//
//
//	delete [] mZ;
//	delete [] mT;
//	delete [] mQ;
//	delete [] mA1;
//	delete [] mP1;
//	delete [] mAHat;
//	delete [] mVHat;
//	delete [] mA;
//	delete [] mP;
//	delete [] mFInv;
//	delete [] mV;
//	delete [] mL;
//}
//
//double LogGamma(double dX, unsigned int iNu)
//{
//	return exp(-dX*iNu-exp(-dX))/tgamma (iNu);
//
//}
//
//double MixtureValue(double dX,  int* iNumOfComp, double* vWeights ,double* vMeans ,double* vVariances  )
//{
//	double dVal=0;
//
//	for(int i=0; i<iNumOfComp[0]; i++)
//	{
//		dVal=dVal+vWeights[i]*gsl_ran_gaussian_pdf (dX-vMeans[i], sqrt(vVariances[i]));
//
//	}
//
//	return dVal;
//}
//
//void TestMixture( int iNu)
//{
//	int iNumOfComp;
//	double* vWeights =new double [10];
//	double* vMeans = new double [10];
//	double* vVariances= new double [10];
//
//
//
//	MixtureLogGamma(  iNu, &iNumOfComp,  vWeights , vMeans ,vVariances);
//
//	double dSum=0;
//	for(int i=0 ; i<iNumOfComp; i++)
//	{
//		cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
//		dSum=dSum+vWeights[i];
//	}
//
//	for(int i=0 ; i<iNumOfComp; i++)
//	{
//		cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//	}
//
//	for(int i=0 ; i<iNumOfComp; i++)
//	{
//		cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//	}
//	cout<<"Sum Weights "<< dSum <<endl;
//
//	double mu = - digammaOnHost(iNu);
//	double sigma = sqrt(trigammaOnHost(iNu));
//
//	double dLow= mu-4*sigma;
//	double dHigh=mu+6*sigma;
//
//	unsigned int iGridSize=50000;
//	double* mGrid=new double[3*(iGridSize+1)];
//
//
//	cout<< "mu "<< mu<<endl;
//	cout<<"sigma "<<sigma<<endl;
//
//	double dError=0;
//	for(int i=0; i<iGridSize+1;i++)
//	{
//		mGrid[3*i]=dLow+i*(dHigh-dLow)/iGridSize;
//		mGrid[3*i+1]=MixtureValue(dLow+i*(dHigh-dLow)/iGridSize,  &iNumOfComp , vWeights , vMeans , vVariances  );
//		mGrid[3*i+2]=LogGamma(dLow+i*(dHigh-dLow)/iGridSize, iNu);
//		dError=dError+fabs(mGrid[3*i+2]-mGrid[3*i+1])/iGridSize;
//		//cout<< fabs(mGrid[3*i+2]-mGrid[3*i+1])<< endl;
//	}
//
//	cout<<"dError "<<dError<<endl;
//
//	string sGridFile="mGrid.csv";
//	WriteOutDoubleArray(mGrid, iGridSize+1, 3, sGridFile);
//	delete [] mGrid;
//	delete [] vWeights;
//	delete [] vMeans;
//	delete [] vVariances;
//
//
//}


void InitializeQuantilePPAlgo(int iN, double dP, double * mMarker ,  double * mN , double * mND)
{
	for(int i=0; i<iN; i++)
	{
		gsl_sort (&mMarker[i*5], 1, 5);

		for(int j=0; j<5; j++)
		{
			mN[i*5+j]=j+1;

		}

		mND[i*5]=1;
		mND[i*5+1]=1+2*dP;
		mND[i*5+2]=1+4*dP;
		mND[i*5+3]=3+2*dP;
		mND[i*5+4]=5;


	}
}

double Sign(double x)
{
	if (x > 0) return 1.0;
	if (x < 0) return -1.0;
	return 0;
}

void QuantilePPAlgo(int iN,   double dP,double* vX, double * mMarker , double *mN , double * mND)
{
	for(int i=0; i<iN;i++)
	{
		/* Find k */

		int k=0;
		if(vX[i]<mMarker[i*5])
		{
			mMarker[i*5]=vX[i];
			k=0;

		}
		else if(vX[i]>mMarker[i*5+4])
		{
			mMarker[i*5+4]=vX[i];
			k=3;

		}
		else
		{

			while( (mMarker[i*5+k+1]< vX[i]))
			{

				k+=1;

			}
		}


		/* Increment N*/
		for(int j=k+1; j<5; j++)
		{
			mN[i*5+j]+=1;
		}
		/* Increment desired N*/
		mND[i*5]+=0;
		mND[i*5+1]+=dP/2;
		mND[i*5+2]+=dP;
		mND[i*5+3]+=(1+dP)/2;
		mND[i*5+4]+=1;

		double  dD;
		/* Adjusting the heights of the markers if necessary */
		for(int j=1; j<4;j++)
		{
			dD=mND[i*5+j]-mN[i*5+j];
			if( ( (dD>=1) & (mN[i*5+j+1]-mN[i*5+j] >1) ) | ( (dD<=-1) & (mN[i*5+j-1]-mN[i*5+j] <-1) ))
			{

				dD=Sign(dD);
				/*  Use P^2 update */
				double dQ=mMarker[i*5+j]+(dD/(mN[i*5+j+1]-mN[i*5+j-1]))*( (mN[i*5+j]-mN[i*5+j-1]+dD)*(mMarker[i*5+j+1]-mMarker[i*5+j])/(mN[i*5+j+1]-mN[i*5+j]) +
						(mN[i*5+j+1]-mN[i*5+j]-dD)*(mMarker[i*5+j]-mMarker[i*5+j-1])/(mN[i*5+j]-mN[i*5+j-1]));
				/* If not okay then use linear update */
				if((dQ>mMarker[i*5+j-1])&(dQ<mMarker[i*5+j+1]))
				{
					mMarker[i*5+j]=dQ;
				}
				else
				{
					if(dD<0)
					{
						mMarker[i*5+j]=mMarker[i*5+j]+dD*(mMarker[i*5+j-1]-mMarker[i*5+j])/(mN[i*5+j-1]-mN[i*5+j]);
					}
					else
					{
						mMarker[i*5+j]=mMarker[i*5+j]+dD*(mMarker[i*5+j+1]-mMarker[i*5+j])/(mN[i*5+j+1]-mN[i*5+j]);
					}
				}

				/* Update N*/
				mN[i*5+j]+=dD;
			}
		}


	}

} 
 

#undef PI

#endif /* BASIC_FUNCTIONS_H_ */
