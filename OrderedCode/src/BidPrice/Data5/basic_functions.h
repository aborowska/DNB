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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <iomanip>

#include "rtnorm.h"


#define PI 3.141592653589793238462643383279502884197169399375105820974944


int isthisnan(double x){return x!=x;}
int isthisinf(double x){return (!isthisnan(x) && isthisnan(x-x)) | (fabs(x)>DBL_MAX );}
int isthisfinite(double x){return !(isthisnan(x) || isthisinf(x));}


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
			vInitialLogVol[i]=log(dVar);
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

double LogitTransform(double x )
{
	return log(x/(1-x));
}

double LogitTransformBack(double x )
{
	return exp(x)/(1+exp(x));
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

void DrawAuxYandAuxHNormal(struct AllParam sParam, int iNumOfObs , double* vData, const gsl_rng * gsl_random_num)
{

	double* vWeights =new double [10];
	double* vMeans = new double [10];
	double* vVariances= new double [10];
	double* dU=new double;
	double* dCumP=new double[10];
	double* dLogVar=new double;
//	double dOffSet=0.000000001;
	double dOffSet=0.000001;

	for(int k=0; k<iNumOfObs; k++)
	{
		if(isthisnan(sParam.vYStar[k])==0 )
		{
	//		cout<<"XXXXXXXXXXXXXXXX  Start XXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;
	//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
	//		cout << "sParam.vN[i] " << sParam.vN[k]<< endl;
	//		cout<< "sParam.vX[i] "<< sParam.vX[k]<<endl;
	//		cout << "sParam.vTau1[k] " << sParam.vTau1[k]<<endl;
	//		cout << "sParam.vTau2[k] " << sParam.vTau2[k]<<endl;

			dLogVar[0]=sParam.dMu[0]+sParam.vS[k]+sParam.vX[k];


			/* Weights */
			vWeights[0]=0.00609; vWeights[1]=0.04775; vWeights[2]=0.13057; vWeights[3]=0.20674; vWeights[4]=0.22715;
			vWeights[5]=0.18842; vWeights[6]=0.12047; vWeights[7]=0.05591; vWeights[8]=0.01575; vWeights[9]=0.00115;

			/* Means */
			vMeans[0]=1.92677; vMeans[1]=1.34744; vMeans[2]= 0.73504; vMeans[3]= 0.02266; vMeans[4]= -0.85173;
			vMeans[5]=-1.97278; vMeans[6]= -3.46788; vMeans[7]=-5.55246; vMeans[8]= -8.68384; vMeans[9]= -14.65000;

			/* Variances */
			vVariances[0]=0.11265; vVariances[1]=0.17788; vVariances[2]=0.26768; vVariances[3]=0.40611; vVariances[4]= 0.62699;
			vVariances[5]=0.98583; vVariances[6]=1.57469; vVariances[7]= 2.54498; vVariances[8]= 4.16591; vVariances[9]= 7.33342;

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

			dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(log(sParam.vYStar[k]*sParam.vYStar[k]+dOffSet)-dLogVar[0]-vMeans[0],2) /vVariances[0]);


			for(int i=1; i<10;i++)
			{
				dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(log(sParam.vYStar[k]*sParam.vYStar[k]+dOffSet)-dLogVar[0]-vMeans[i],2) /vVariances[i]);
			}




			int i=0;
			dU[0]=gsl_rng_uniform(gsl_random_num);

			while((dCumP[i]/dCumP[10-1])<dU[0])
			{
				i=i+1;

			}


			/*Saving AuxY and AuxH for Tau1  */



			sParam.mAuxY[k]=log(sParam.vYStar[k]*sParam.vYStar[k]+dOffSet)-vMeans[i];
			sParam.mAuxH[k]=vVariances[i];

	//		cout<<"Mixture component at " << k << " is "<< i<<endl;
	//		cout<<"AuxY at " << k << " is "<< sParam.mAuxY[k]<<endl;
	//		cout<<"AuxH at " << k << " is "<< sParam.mAuxH[k]<<endl;

		}
		else
		{
			sParam.mAuxY[k]=sParam.vYStar[k];
			sParam.mAuxH[k]=sParam.vYStar[k];
		}
//		cout<<"AuxY at " << k << " is "<< sParam.mAuxY[k]<<endl;
//		cout<<"YStart at " << k << " is "<< sParam.vYStar[k]<<endl;
	}

	delete [] vWeights;
	delete [] vMeans;
	delete [] vVariances;
	delete [] dCumP;

	delete dU;
	delete dLogVar;

}

void DrawAuxYandAuxHT(struct AllParam sParam, int iNumOfObs , double* vData, const gsl_rng * gsl_random_num)
{

	double* vWeights =new double [10];
	double* vMeans = new double [10];
	double* vVariances= new double [10];
	double* dU=new double;
	double* dCumP=new double[10];
	double* dLogVar=new double;
//	double dOffSet=0.000000001;
	double dOffSet=0.001;

	for(int k=0; k<iNumOfObs; k++)
	{
//		if(isthisnan(sParam.vYStar[k])==0 )
//		{
	//		cout<<"XXXXXXXXXXXXXXXX  Start XXXXXXXXXXXXXXXXXXXXXXXxx "<< endl;
	//		cout << " sParam.dMu[0] "<<  sParam.dMu[0]<< endl;
	//		cout << "sParam.vN[i] " << sParam.vN[k]<< endl;
	//		cout<< "sParam.vX[i] "<< sParam.vX[k]<<endl;
	//		cout << "sParam.vTau1[k] " << sParam.vTau1[k]<<endl;
	//		cout << "sParam.vTau2[k] " << sParam.vTau2[k]<<endl;

			dLogVar[0]=sParam.dMu[0]+sParam.vS[k]+sParam.vX[k];


			/* Weights */
			vWeights[0]=0.00609; vWeights[1]=0.04775; vWeights[2]=0.13057; vWeights[3]=0.20674; vWeights[4]=0.22715;
			vWeights[5]=0.18842; vWeights[6]=0.12047; vWeights[7]=0.05591; vWeights[8]=0.01575; vWeights[9]=0.00115;

			/* Means */
			vMeans[0]=1.92677; vMeans[1]=1.34744; vMeans[2]= 0.73504; vMeans[3]= 0.02266; vMeans[4]= -0.85173;
			vMeans[5]=-1.97278; vMeans[6]= -3.46788; vMeans[7]=-5.55246; vMeans[8]= -8.68384; vMeans[9]= -14.65000;

			/* Variances */
			vVariances[0]=0.11265; vVariances[1]=0.17788; vVariances[2]=0.26768; vVariances[3]=0.40611; vVariances[4]= 0.62699;
			vVariances[5]=0.98583; vVariances[6]=1.57469; vVariances[7]= 2.54498; vVariances[8]= 4.16591; vVariances[9]= 7.33342;

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

			dCumP[0]=exp(log(vWeights[0])-0.5*log(vVariances[0])-0.5*pow(log(sParam.vYStar[k]*sParam.vYStar[k]+dOffSet)-log(sParam.vLambda[k])-dLogVar[0]-vMeans[0],2) /vVariances[0]);


			for(int i=1; i<10;i++)
			{
				dCumP[i]=dCumP[i-1]+exp(log(vWeights[i])-0.5*log(vVariances[i])-0.5*pow(log(sParam.vYStar[k]*sParam.vYStar[k]+dOffSet)-log(sParam.vLambda[k])-dLogVar[0]-vMeans[i],2) /vVariances[i]);
			}




			int i=0;
			dU[0]=gsl_rng_uniform(gsl_random_num);

			while((dCumP[i]/dCumP[10-1])<dU[0])
			{
				i=i+1;

			}


			/*Saving AuxY and AuxH for Tau1  */



			sParam.mAuxY[k]=log(sParam.vYStar[k]*sParam.vYStar[k]+dOffSet)-log(sParam.vLambda[k])-vMeans[i];
			sParam.mAuxH[k]=vVariances[i];

	//		cout<<"Mixture component at " << k << " is "<< i<<endl;
	//		cout<<"AuxY at " << k << " is "<< sParam.mAuxY[k]<<endl;
	//		cout<<"AuxH at " << k << " is "<< sParam.mAuxH[k]<<endl;

//		}
//		else
//		{
//			sParam.mAuxY[k]=sParam.vYStar[k];
//			sParam.mAuxH[k]=sParam.vYStar[k];
//		}
//		cout<<"AuxY at " << k << " is "<< sParam.mAuxY[k]<<endl;
//		cout<<"YStart at " << k << " is "<< sParam.vYStar[k]<<endl;
	}

	delete [] vWeights;
	delete [] vMeans;
	delete [] vVariances;
	delete [] dCumP;

	delete dU;
	delete dLogVar;

}

void GammaLogPosterior(double * dLogGamma, double *dPriorA, double *dPriorB,  int iNumOfObs,
		double * vProb, double* vIndicator, double * dLogPosterior)
{
	dLogPosterior[0]=0;
	double dGamma=LogitTransformBack(dLogGamma[0]);

	for(int i=0; i<iNumOfObs; i++)
	{
		dLogPosterior[0]+=log( dGamma*vIndicator[i]+(1-dGamma)*vProb[i] );
	}
	dLogPosterior[0]+=(dPriorA[0] -1)*log(dGamma)+(dPriorB[0] -1)*log(1-dGamma);
}

void DrawGammaAdaptiveRWNormal(  double* vData, int iNumOfObs,  int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam,
		 double* dCovar,  double* dSum  )
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


	double * vNormal =new double[iNumOfObs];
	double * vIndicator = new double[iNumOfObs];
	for(int i=0; i<iNumOfObs; i++)
	{
		if(vData[i]==0)
		{
			vIndicator[i]=1;
		}
		else
		{
			vIndicator[i]=0;
		}

		vNormal[i]=gsl_cdf_gaussian_P(vData[i]+0.5, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) )-gsl_cdf_gaussian_P(vData[i]-0.5, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) );

	}

	GammaLogPosterior(&dNewGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vNormal, vIndicator, dLogPosteriorNew)	;
	GammaLogPosterior(&dOldGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vNormal, vIndicator, dLogPosteriorOld)	;

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
		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));
	}

//	cout << "dSum " << dSum[0]<< endl;
//	cout << "dCovar " << dCovar[0]<< endl;
//	cout << "proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
//	cout << "accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;
	delete [] vNormal;
	delete [] vIndicator;

}

void DrawGammaAdaptiveRWT(  double* vData, int iNumOfObs,  int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam,
		 double* dCovar,  double* dSum  )
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


	double * vNormal =new double[iNumOfObs];
	double * vIndicator = new double[iNumOfObs];
	for(int i=0; i<iNumOfObs; i++)
	{
		if(vData[i]==0)
		{
			vIndicator[i]=1;
		}
		else
		{
			vIndicator[i]=0;
		}

		vNormal[i]=gsl_cdf_tdist_P((vData[i]+0.5)/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), sParam.dNu[0] )-gsl_cdf_tdist_P((vData[i]-0.5)/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), sParam.dNu[0] );

	}

	GammaLogPosterior(&dNewGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vNormal, vIndicator, dLogPosteriorNew)	;
	GammaLogPosterior(&dOldGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vNormal, vIndicator, dLogPosteriorOld)	;

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
		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1));
	}

//	cout << "dSum " << dSum[0]<< endl;
//	cout << "dCovar " << dCovar[0]<< endl;
//	cout << "proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
//	cout << "accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;
	delete [] vNormal;
	delete [] vIndicator;

}

double DrawTruncatedNormal(double dMu, double  dSigma, double  dLow, double dHigh ,const gsl_rng * gsl_random_num)
{
	double dLowS=dLow/dSigma;
	double dHighS=dHigh/dSigma;


	int iExit=0;
	double dDraw,dPhi;

	while(iExit==0)
	{
		double dZ=gsl_rng_uniform(gsl_random_num)*(dHighS-dLowS)+dLowS;

		if((dLowS<0)&(dHighS>0))
		{
			dPhi=exp(-dZ*dZ/2);
		}
		else if(dLowS >0)
		{
			dPhi=exp((dLowS*dLowS-dZ*dZ)/2);
		}
		else
		{
			dPhi=exp((dHighS*dHighS-dZ*dZ)/2);
		}

		double dU=gsl_rng_uniform(gsl_random_num);
		if(dU<=dPhi)
		{
			dDraw=dZ*dSigma;
			iExit=1;
		}
	}

	return dDraw;
}

double DrawTruncatedNormalOLD(double dMu, double  dSigma, double  dLow, double dHigh ,const  gsl_rng * gsl_random_num)
{
//	const gsl_rng_type * T=gsl_rng_default;
//	gsl_rng *random_num=gsl_rng_alloc(T);


	std::pair<double, double> s;

//	const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
//	gsl_rng *gen = gsl_rng_alloc (type);
//	s = rtnorm(gsl_random_num,dLow,dHigh,dMu,dSigma);
	return s.first;
}

double CalcVariance( int iNumOfObs, double * vX)
{
	double dMeanX=0;
	double dMeanX2=0;

	for(int i=0; i<iNumOfObs; i++)
	{
		dMeanX+=vX[i]/iNumOfObs;
		dMeanX2+=(vX[i]*vX[i])/iNumOfObs;
	}

	return dMeanX2-dMeanX*dMeanX;
}

double DrawTruncatedTOLD( double  dSigma, double dNu, double  dLow, double dHigh , const gsl_rng * gsl_random_num)
{
//	double dDraw;
//	double dU=gsl_ran_chisq(gsl_random_num, dNu);
//
//	double dScale=sqrt(dU/dNu)/dSigma;
//	dDraw=DrawTruncatedNormal(0, 1, dLow*dScale, dHigh*dScale , gsl_random_num);
//
//	return dDraw/dScale;

	double dDraw;
	double dScale=dSigma*sqrt(1/gsl_ran_gamma(gsl_random_num, dNu/2 , 2/dNu));
	dDraw=DrawTruncatedNormal(0, 1, dLow/dScale, dHigh/dScale , gsl_random_num);

	return dDraw*dScale;


}

double DrawTruncatedT( double  dSigma, double dNu, double  dLow, double dHigh , const gsl_rng * gsl_random_num)
{
	double dU=gsl_rng_uniform(gsl_random_num);

	return gsl_cdf_tdist_Pinv(gsl_cdf_tdist_P(dLow/dSigma,dNu)+dU*(gsl_cdf_tdist_P(dHigh/dSigma,dNu)-gsl_cdf_tdist_P(dLow/dSigma,dNu) ) , dNu)*dSigma;


}


void DrawYStarNormal(struct AllParam sParam,  double* vData,int iNumOfObs , const gsl_rng * gsl_random_num)
{
	double dU;
	double dThresh;
//	pair<double, double> s;
//	 gsl_rng_env_setup();                          // Read variable environnement
//	  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
//	  gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation

	for(int i=0; i<iNumOfObs; i++)
	{
		if(vData[i]==0)
		{
			dU=gsl_rng_uniform(gsl_random_num);
			dThresh=sParam.dGamma[0]/(sParam.dGamma[0]+(1-sParam.dGamma[0])*(gsl_cdf_gaussian_P(0.5, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) )-gsl_cdf_gaussian_P(-0.5, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) ) ) );

			if(dU<=dThresh)
			{
//				cout <<"dU "<<dU<<endl;
//				cout <<" Gamma "<< sParam.dGamma[0]<<endl;


				sParam.vYStar[i]=gsl_ran_gaussian(gsl_random_num, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2));

			}
			else
			{
				sParam.vYStar[i]=DrawTruncatedNormal(0, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), -0.5, 0.5 , gsl_random_num);

//				s = rtnorm( gen,-0.5, 0.5 ,0,exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2));
//				sParam.vYStar[i]=s.first;
//				cout<<"sParam.vYStar[i]"<< sParam.vYStar[i]<<endl;
			}
//				sParam.vYStar[i]=vYStarTrue[i];
		}
		else
		{
			sParam.vYStar[i]=DrawTruncatedNormal(0, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), vData[i]-0.5, vData[i]+0.5 , gsl_random_num);

//			s = rtnorm( gen,vData[i]-0.5,vData[i]+ 0.5 ,0,exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2));
//			sParam.vYStar[i]=s.first;
		}
//		cout <<"Data and vY Star at "<<i<<" is "<<vData[i] <<" and  "<< sParam.vYStar[i]<< endl;
	}

//	cout << "Var y*  is " <<CalcVariance( iNumOfObs, sParam.vYStar)<<endl;
}


void DrawYStarT(struct AllParam sParam,  double* vData,int iNumOfObs , const gsl_rng * gsl_random_num)
{
	double dU;
	double dThresh;

	for(int i=0; i<iNumOfObs; i++)
	{
		if(vData[i]==0)
		{
			dU=gsl_rng_uniform(gsl_random_num);
			dThresh=sParam.dGamma[0]/(sParam.dGamma[0]+(1-sParam.dGamma[0])*(gsl_cdf_gaussian_P(0.5/sqrt(sParam.vLambda[i]), exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) )-gsl_cdf_gaussian_P(-0.5/sqrt(sParam.vLambda[i]), exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) ) ) );

			if(dU<=dThresh)
			{
//				cout <<"dU "<<dU<<endl;
//				cout <<" Gamma "<< sParam.dGamma[0]<<endl;


				sParam.vYStar[i]=sqrt(sParam.vLambda[i])*gsl_ran_gaussian(gsl_random_num, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2));

			}
			else
			{
				sParam.vYStar[i]=sqrt(sParam.vLambda[i])*DrawTruncatedNormal(0, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), -0.5/sqrt(sParam.vLambda[i]), 0.5/sqrt(sParam.vLambda[i]) , gsl_random_num);


			}

		}
		else
		{
			sParam.vYStar[i]=sqrt(sParam.vLambda[i])*DrawTruncatedNormal(0, exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), (vData[i]-0.5)/sqrt(sParam.vLambda[i]), (vData[i]+0.5)/sqrt(sParam.vLambda[i]) , gsl_random_num);


		}

	}


}

void DrawYStarTOLD(struct AllParam sParam,  double* vData,int iNumOfObs , const gsl_rng * gsl_random_num)
{
	double dU;
	double dThresh;

//	double dRoundT=0.000000005;
	double dRoundT=0.5;
	for(int i=0; i<iNumOfObs; i++)
	{
		if(vData[i]==0)
		{
			dU=gsl_rng_uniform(gsl_random_num);
			dThresh=sParam.dGamma[0]/(sParam.dGamma[0]+(1-sParam.dGamma[0])*(gsl_cdf_tdist_P(dRoundT/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), sParam.dNu[0]  )
										-gsl_cdf_tdist_P(-dRoundT/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2), sParam.dNu[0]  ) ) );

			if(dU<=dThresh)
			{
//				cout <<"dU "<<dU<<endl;
//				cout <<" Gamma "<< sParam.dGamma[0]<<endl;


				sParam.vYStar[i]=gsl_ran_tdist(gsl_random_num, sParam.dNu[0])*exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2);//*sqrt((sParam.dNu[0]-2)/sParam.dNu[0]);

			}
			else
			{
				sParam.vYStar[i]=DrawTruncatedT( exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2),sParam.dNu[0], -dRoundT, dRoundT , gsl_random_num);



			}

		}
		else
		{
			sParam.vYStar[i]=DrawTruncatedT( exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2),sParam.dNu[0], vData[i]-dRoundT, vData[i]+dRoundT , gsl_random_num);


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

	gsl_matrix_free (mPUnique);
	gsl_matrix_free (mWUnique);
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
	gsl_matrix_free (mLambdaInvTheta);

	gsl_permutation_free (mPermutation);

	delete dLambda;
	delete dH;
	delete dHPlus;

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

void SystemMatricesNormalWithSeasonal( int iNumOfKnots, struct AllParam sParam,   double* mZ, double* mT, double *mQ, double* mA1, double* mP1)
{
	/* Z */
	for(int i=0; i<iNumOfKnots-1;i++)
	{
		mZ[i]=0;

	}
	for(int i=0; i<2;i++)
	{
		mZ[iNumOfKnots-1+i]=1;

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

void KalmanFilterNormalWithSeasonal(int iNumOfKnots, double dGamma,  gsl_matrix * mWTilde , double* mY,double * vData, double* mH, double* mZ, double* mT, double *mQ,
		double* mA1, double* mP1, int iDimOfObs, int iDimOfStates,  int iNumOfObs, double* mA, double * mP, double* mFInv, double* mV, double* mL)
{

	double * mZTrans= new double[iDimOfObs*iDimOfStates];
	double * mKv= new double[iDimOfStates];
	double * mKZ= new double[iDimOfStates*iDimOfStates];
	double * mTa= new double[iDimOfStates];
	double * mLTrans=new double[iDimOfStates*iDimOfStates];
	double * mPL=new double[iDimOfStates*iDimOfStates];
	double * mTPL=new double[iDimOfStates*iDimOfStates];

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

	double * mZa= new double;
	double * mZPZ= new double;
	double * mK= new double[iDimOfObs*iDimOfStates];
	double * mPZFInv = new double[iDimOfStates*iDimOfObs];
	double * mZFInv = new double[iDimOfStates*iDimOfObs];
	double * mF= new double;






	for(int i=0; i<iNumOfObs; i++)
	{
		/* Seting the time varying Z matrix */
		for(int k=0; k<iNumOfKnots-1;k++)
		{
				mZ[k]=gsl_matrix_get(mWTilde,i,k);

		}
//		for(int l=0; l<(iNumOfKnots-1+2)*2; l++)
//		{
//				cout<< "At obs "<<i <<" mZ at "<<l<<" is "<<mZ[l]<<endl;
//		}




		MatrixTrans( mZ, iDimOfObs, iDimOfStates, mZTrans);

//		for(int l=0; l<iDimOfStates*iDimOfObs; l++)
//		{
//			cout<< "At obs "<<i <<" mZ at "<<l<<" is "<<mZTrans[l]<<endl;
//		}



		/* v_t */
//		if(isthisnan(vData[i])!=0 )
//		{
//			mV[0]=0;
////			cout<<"nan y  at "<< i << " is "<<mY[i]<<endl;
////			cout<<"nan data "<< i << " is "<<vData[i]<<endl;
//		}
//		else
//		{
			if(i==0)
			{
				MatrixMulti(mZ, mA1,iDimOfObs,iDimOfStates, iDimOfStates, 1, mZa);
			}
			else
			{
				MatrixMulti(mZ, mA,iDimOfObs,iDimOfStates, iDimOfStates, 1, mZa);
			}
			MatrixSub(mY, mZa, iDimOfObs, 1, mV);
//		}
//		for(int l=0; l<iDimOfObs; l++)
//		{
//			cout<< "At obs "<<i <<" mZa at "<<l<<" is "<<mZa[l]<<endl;
//		}
//
//		for(int l=0; l<iDimOfObs; l++)
//		{
//			cout<< "At obs "<<i <<" mY at "<<l<<" is "<<mY[l]<<endl;
//		}



//		for(int l=0; l<iDimOfObs; l++)
//		{
//			cout<< "At obs "<<i <<" mV at "<<l<<" is "<<mV[l]<<endl;
//		}


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
//		if(isthisnan(vData[i])!=0 )
//		{
//			for(int l=0; l<iDimOfStates;l++)
//			{
//				for(int m=0; m<iDimOfObs;m++)
//				{
//					mK[l*iDimOfObs+m]=0;
//				}
//			}
//			mFInv[0]=0;
//		}
//		else
//		{
			MatirxInv2x2(mF , iDimOfObs,mFInv);
			MatrixMulti(mZTrans, mFInv,iDimOfStates, iDimOfObs, iDimOfObs, iDimOfObs, mZFInv);
			MatrixMulti(mP, mZFInv,iDimOfStates, iDimOfStates, iDimOfStates, iDimOfObs, mPZFInv);
			MatrixMulti(mT, mPZFInv,iDimOfStates, iDimOfStates, iDimOfStates, iDimOfObs, mK);
//		}
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





	delete mZa;
	delete mZPZ;
	delete  mF;
	delete [] mKZ;
	delete [] mK;
	delete [] mPZFInv;
	delete [] mZFInv;
	delete [] mZTrans;
	delete [] mTa;
	delete [] mKv;
	delete [] mPL;
	delete [] mTPL;
	delete [] mLTrans;


}

void KalmanSmootherNormalWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde , double* mZ,  double * mA, double *mP, double *mV,
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
		}



		MatrixTrans( mZ, iDimOfObs, iDimOfStates, mZTrans);
		/* mR */
		MatrixTrans( mL, iDimOfStates, iDimOfStates, mLTrans);
//		for(int k=0; k<iDimOfStates; k++)
//		{
//			for(int j=0; j<iDimOfStates; j++)
//			{
//				cout<<"mL"<< k << ","<<j <<" at "<<i<<" is "<<mL[k*iDimOfStates+j]<<endl;
//			}
//
//		}
		MatrixMulti(mFInv, mV,iDimOfObs, iDimOfObs, iDimOfObs, 1, mFv);

//		cout<<"mFInv "<<mFInv[0]<<endl;
		MatrixMulti(mZTrans, mFv,iDimOfStates, iDimOfObs, iDimOfObs, 1, mZFv);
		MatrixMulti(mLTrans, mR,iDimOfStates, iDimOfStates, iDimOfStates,1, mLR);

//		for(int k=0; k<iDimOfStates; k++)
//		{
//			cout<<"mZFv at "<<i<<" is "<< mZFv[k]<<endl;
//		}

		MatrixAdd( mZFv, mLR, iDimOfStates, 1, mNextR);

//		for(int k=0; k<iDimOfStates; k++)
//		{
//
//			cout<<"R at "<<i<<" is "<< mNextR[k]<<endl;
//		}

		/* AHat */
		MatrixMulti(mP, mNextR,iDimOfStates, iDimOfStates, iDimOfStates,1, mPR);
		MatrixAdd( mA, mPR, iDimOfStates, 1, mAHat);

//		for(int k=0; k<iDimOfStates; k++)
//		{
//			for(int j=0; j<iDimOfStates; j++)
//			{
//						cout<<"P"<< k << ","<<j <<" at "<<i<<" is "<<mP[k*iDimOfStates+j]<<endl;
//			}
//
//		}
//
//		for(int k=0; k<iDimOfStates; k++)
//		{
//			cout<<"AHat at "<<i<<" is "<<mAHat[k]<<endl;
//		}

		/* mN*/
		MatrixSandwitch(mZ, mFInv,iDimOfObs,iDimOfStates, iDimOfObs, iDimOfObs, mZFZ);
		MatrixSandwitch(mL, mN,iDimOfStates,iDimOfStates, iDimOfStates, iDimOfStates, mLNL);
		MatrixAdd( mZFZ, mLNL, iDimOfStates, iDimOfStates, mNextN);

		/* VHat */
		MatrixSandwitch(mP, mNextN,iDimOfStates,iDimOfStates, iDimOfStates, iDimOfStates, mPNP);
		MatrixSub( mP, mPNP, iDimOfStates, iDimOfStates, mVHat);




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

void CalculateLLNormalWithSeasonal( int iNumOfKnots, double* vData,  gsl_matrix * mWTilde ,struct AllParam sParam,int iNumOfObs,  int iDimOfObs,
		int iDimOfStates, double * dLL )
{





	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatricesNormalWithSeasonal(iNumOfKnots, sParam, mZ, mT, mQ,  mA1,  mP1);

	double * mA =new double[iNumOfObs*iDimOfStates];
	double * mP =new double[iNumOfObs*iDimOfStates*iDimOfStates];
	double * mFInv =new double[iNumOfObs*iDimOfObs*iDimOfObs];
	double * mV =new double[iNumOfObs*iDimOfObs];
	double * mL=new double[iNumOfObs*iDimOfStates*iDimOfStates];
	KalmanFilterNormalWithSeasonal(iNumOfKnots, sParam.dGamma[0],  mWTilde ,sParam.mAuxY, vData, sParam.mAuxH,   mZ, mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);

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


		Determinant2x2(mFInv,iDimOfObs, dDetFInv);
		MatrixSandwitch(mV, mFInv,iDimOfObs,1, iDimOfObs, iDimOfObs, dvFv);
		dLL[0]=dLL[0]-0.5*iDimOfObs*log(2*PI)+0.5*(log(dDetFInv[0])-dvFv[0]);

//		cout << "LL Cont at " << i << " is " << 0.5*(log(dDetFInv[0])-dvFv[0]) << endl;
//		cout << "LL at " << i << " is " << dLL[0] << endl;


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

void SSFRecursionNormalWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam, double *mH, double* mZ, double* mT, double *mQ, double* mA1, double* mP1,
		int iNumOfObs,  int iDimOfObs, int iDimOfStates, const gsl_rng * gsl_random_num, double* mYPlus, double* mAPlus)
{



	double * mQChol =new double[iDimOfStates*iDimOfStates];
//	double * mHChol =new double[iDimOfObs*iDimOfObs];
	double * mHChol =new double;
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


		}


//		cout << "SSF Recuriosn H Chol START " <<endl;
//		cout<<"vN[i] "<<vN[i]<<endl;

//		Cholesky(mH, iDimOfObs, mHChol);

		mHChol[0]=sqrt(mH[0]);
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

//		cout<<"Recursion mEps"<<mEps[0]<<endl;
//		cout<<"Recursion mZa"<<mZa[0]<<endl;
//		cout << "mAPlus[0] " <<mNextAPlus[0] <<endl;
//		cout << "mAPlus[1] " <<mNextAPlus[1] <<endl;
//		cout << "mAPlus[2] " <<mNextAPlus[2] <<endl;
//		cout << "mAPlus[3] " <<mNextAPlus[3] <<endl;
//		cout<<"Recursion Y"<<mNextYPlus[0]<<endl;

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
//
//		}
//		else
//		{
//			for(int k=0; k<2;k++)
//			{
//				mRandom[(2+iDimOfStates)*(i+1)+iDimOfStates+k]=mEps[k];
//			}
//
//		}




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




	delete [] mQChol ;
	delete  mHChol;
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

void SimulationSmootherNormalWithSeasonal( int iNumOfKnots, double* vData, gsl_matrix * mWTilde ,struct AllParam sParam,double *mY,double *mH,unsigned int iNumOfObs,
		unsigned int iDimOfObs,unsigned int iDimOfStates, const gsl_rng * gsl_random_num, double* mDraw)
{
	double * mZ =new double[iDimOfObs*iDimOfStates];
	double * mT =new double[iDimOfStates*iDimOfStates];
	double * mQ =new double[iDimOfStates*iDimOfStates];
	double * mA1 =new double[iDimOfStates];
	double * mP1 =new double[iDimOfStates*iDimOfStates];
	SystemMatricesNormalWithSeasonal(iNumOfKnots,sParam,mZ, mT, mQ, mA1, mP1);

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


	SSFRecursionNormalWithSeasonal(iNumOfKnots,    mWTilde ,sParam,mH, mZ,  mT, mQ, mA1,  mP1,  iNumOfObs,  iDimOfObs,iDimOfStates,  gsl_random_num, mYPlus, mAPlus);

//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mAPlus at " << i<< " is "<< mAPlus[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfObs; i++)
//	{
//		cout <<"mYPlus at " << i<< " is "<< mYPlus[i]<< endl;
//	}


	/* AHat */
	KalmanFilterNormalWithSeasonal(iNumOfKnots, sParam.dGamma[0],  mWTilde ,mY, vData, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
	KalmanSmootherNormalWithSeasonal( iNumOfKnots,    mWTilde ,mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHat,  mVHat );

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
	KalmanFilterNormalWithSeasonal(iNumOfKnots,  sParam.dGamma[0],   mWTilde ,mYPlus, vData, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mA at " << i<< " is "<< mA[i]<< endl;
//	}
//
//	for(int i=0 ; i<iNumOfObs*iDimOfStates; i++)
//	{
//		cout <<"mV at " << i<< " is "<< mV[i]<< endl;
//	}

	KalmanSmootherNormalWithSeasonal( iNumOfKnots,    mWTilde , mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHatPlus,  mVHatPlus);

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

void DrawXandMuWithSeasonal( int iNumOfKnots,   gsl_matrix * mWTilde ,struct AllParam sParam, int iNumOfObs,const gsl_rng * gsl_random_num )
{
	int iDimOfObs=1;
	int iDimOfStates= iNumOfKnots-1+2;
	double * mDraw=new double[ iDimOfStates*iNumOfObs];

//	for(int i=0;i<iNumOfObs; i++)
//	{
//		cout<<"Draw state sParam.vYStar at "<<i<<" is "<<sParam.vYStar[i]<<endl;
//		cout<<"Draw state sParam.mAuxY at "<<i<<" is "<<sParam.mAuxY[i]<<endl;
//	}


	SimulationSmootherNormalWithSeasonal(iNumOfKnots,sParam.vYStar, mWTilde ,sParam,sParam.mAuxY,sParam.mAuxH, iNumOfObs,iDimOfObs,iDimOfStates, gsl_random_num, mDraw);


	for(int j=0; j<iNumOfKnots-1; j++)
	{
		cout <<"beta " <<mDraw[j]<<endl;
	}

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

void PhiSigmaLogPosteriorWithSeasonal(int iNumOfKnots, double* vData,  gsl_matrix * mWTilde ,struct AllParam sParam, int iNumOfObs, int iDimOfObs,
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
	CalculateLLNormalWithSeasonal( iNumOfKnots, vData, mWTilde , sParam,iNumOfObs, iDimOfObs, iDimOfStates,  dLL );

	dLogPost[0]=dLogPost[0]+dLL[0];
//	cout<<"dLL "<<dLL[0]<<endl;
//	cout<<"full posterior "<<dLogPost[0]<<endl;
	delete dLL;

}

void DrawPhiSigmaAdaptiveRWWithSeasonal( int iNumOfKnots, gsl_matrix * mWTilde ,struct AllParam sParam,int iNumOfObs,  int iNumOfIter,
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
	PhiSigmaLogPosteriorWithSeasonal(iNumOfKnots, sParam.vYStar,  mWTilde,sParam, iNumOfObs, 1, iNumOfKnots-1+2,dNewLogPosterior);

	sParam.dPhi[0]=LogitTransformBack(mLogParam[0]);
	sParam.dSigma2[0]=exp(mLogParam[1]);

//	cout<<"XXXXX old XXXXXXXX"<<endl;
//	cout<<" old Phi " << sParam.dPhi[0]<< endl;
//	cout<<" old Sigma2 " << sParam.dSigma2[0]<< endl;
	PhiSigmaLogPosteriorWithSeasonal(iNumOfKnots,  sParam.vYStar, mWTilde,sParam, iNumOfObs, 1, iNumOfKnots-1+2,dOldLogPosterior);

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

void DrawNuAndLambdaT(struct AllParam sParam ,  int iNumOfObs, const gsl_rng * gsl_random_num, double dResolution,  double dDelta)
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

	cout<< "dNuProposed "<<dNuProposed<<endl;
	cout<< "dNuCurrent "<<sParam.dNu[0]<<endl;

	if((dNuProposed < 128) & (dNuProposed >2))
	{
		/* Calculate acceptance probability */
		double dLogProposal=0;
		double dLogCurrent=0;
		for(int i=0; i<iNumOfObs; i++)
		{
			dLogProposal+=log(gsl_ran_tdist_pdf(sParam.vYStar[i]/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2),dNuProposed) );
			dLogCurrent+=log(gsl_ran_tdist_pdf(sParam.vYStar[i]/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i])/2) ,sParam.dNu[0]));
		}

		cout<< "dLogProposed "<<dLogProposal<<endl;
		cout<< "dLogCurrent "<<dLogCurrent<<endl;

		double dAccept=fmin(dLogProposal-dLogCurrent,0);

		/* Accept or reject*/
		double dLogU=log(gsl_rng_uniform(gsl_random_num));

		cout << "dAccept "<<dAccept<<" dLogU "<< dLogU<<endl;
		if(dLogU<=dAccept)
		{
			sParam.dNu[0]=dNuProposed;
		}
	}

	for(int i=0; i<iNumOfObs; i++)
	{

		sParam.vLambda[i]=1.0/gsl_ran_gamma(gsl_random_num, (sParam.dNu[0]+1)/2,2/(sParam.dNu[0]+(sParam.vYStar[i]*sParam.vYStar[i])/exp((sParam.dMu[0]+sParam.vS[i]+sParam.vX[i]))));
	}

}

double OnlineMean(int iN, double dNewObs, double dLastMean )
{
	return (dLastMean+(dNewObs-dLastMean)/iN);
}

double OnlineVar(int iN, double dNewObs, double dLastMean, double dLastVar )
{
	double dNewMean=OnlineMean( iN,  dNewObs,  dLastMean );

	return (dLastVar+dLastMean*dLastMean-dNewMean*dNewMean+(dNewObs*dNewObs-dLastVar-dLastMean*dLastMean)/iN);
}

double VarOrdNormal(double dSigma)
{
	double dT=0.00001;
	double dCont=1;
	double dVar=0;
	double i=1;
	while(dCont>dT)
	{
		dCont=2*i*i*(gsl_cdf_gaussian_P(i+0.5,dSigma) - gsl_cdf_gaussian_P(i-0.5,dSigma));
		dVar+=dCont;
		i+=1;

	}

	return dVar;
}

double VarOrdT(double dSigma, double dNu)
{
	double dT=0.00001;
	double dCont=1;
	double dVar=0;
	double i=1;
	while(dCont>dT)
	{
		dCont=2*i*i*(gsl_cdf_tdist_P((i+0.5)/dSigma, dNu) - gsl_cdf_tdist_P((i-0.5)/dSigma, dNu));
		dVar+=dCont;
		i+=1;

	}

	return dVar;
}

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

void SequentialStore(int iIter, int iMaxIter, int iN, double * vX ,  fstream&   fs ,string sFile)
{
	if(iIter==0)
	{
		fs.open(("OutPut/"+sFile).c_str(), ios::out);
	}
	for(int i=0; i<iN; i++)
	{
		if(i==iN-1)
		{
			fs << vX[i]  << endl;
		}
		else
		{
			fs << vX[i]  << ",";
		}
	}

	if(iIter==iMaxIter-1)
	{
		fs.close();
	}
}

#undef PI

#endif /* BASIC_FUNCTIONS_H_ */
