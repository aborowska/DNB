/*
 * basic_function2.h
 *
 *  Created on: Aug 19, 2014
 *      Author: istvan
 */

#ifndef BASIC_FUNCTION2_H_
#define BASIC_FUNCTION2_H_
#include "struct.h"
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <iomanip>
#include <list>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>

#define PI 3.141592653589793238462643383279502884197169399375105820974944

double HyperGeoOut(double a, double b, double c, double z)
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
//		cout<< i<<endl;
	}
	return S;
}

void CreatPrefixOut(string sInput,string sType, string *sPrefix )
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


void CalculateWStarOut(int iNumOfObs, double* vTimes,  int iNumOfKnots, double* vKnots,gsl_matrix * mLambdaInvTheta,  gsl_matrix * mWStar )
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

void CalculateSplineMatricesOut(int iNumOfObs, double* vTimes,  int iNumOfKnots, double* vKnots, gsl_matrix * mWTilde)
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

	/* WStar matrix */
	cout <<"calculate star started "<<endl;
	CalculateWStarOut( iNumOfObs, vTimes, iNumOfKnots,  vKnots, mLambdaInvTheta, mWStar );
	cout <<"calculate star finished "<<endl;

	/* WTilde matrix */
	for(int i=0; i<iNumOfObs; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			gsl_matrix_set (mWTilde ,i,j, gsl_matrix_get (mW ,i,j) - gsl_matrix_get (mW ,i,iNumOfKnots-1)*gsl_matrix_get (mWStar ,0,j)/gsl_matrix_get (mWStar ,0,iNumOfKnots-1) );
		}
	}

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

void WriteOutDoubleArrayOut(double * mData,int iRows,  int iCols, string  sFile )
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

void ImportMatrix(string sFile,  int iBurnIn,  int iThin, double** mMatrix, double** vMean, int* iNumOfParam, int* iNumOfVar)
{
	vector<vector<double> > mData;

	fstream ifsDataFile;
	ifsDataFile.open(sFile.c_str(), ios::in);
	if(!ifsDataFile.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	string line;
	int iNumOfRows=0;
	int iNumOfColumns=0;
	while ( getline (ifsDataFile,line) )
	{
		vector<double>vRow;
		stringstream ssLine(line);
		string var;
		char* end;
		while(getline(ssLine, var,',') )
		{
			if(iNumOfRows==0)
			{
				iNumOfColumns=iNumOfColumns+1;
			}
			vRow.push_back(strtod( var.c_str() ,&end) );
		}
		mData.push_back(vRow);

		iNumOfRows=iNumOfRows+1;
		vRow.clear();
	}

	ifsDataFile.close();

	int iNumOfMatrixRows=(iNumOfRows-iBurnIn)/iThin;
	int iNextIndex=iBurnIn+iThin;
	int iCurrentMatrixIndex=0;
	*mMatrix = new double[ iNumOfMatrixRows*iNumOfColumns];
	*vMean=new double [iNumOfColumns];

	for(int j=0; j<iNumOfColumns;j++)
	{
		(*vMean)[j]=0;
	}

	for(vector<vector<double> >::iterator i=(mData.begin()+iBurnIn); i != mData.end(); ++i)
	{
		for(vector<double>::iterator j= i->begin(); j != i -> end(); ++j)
		{
			(*vMean)[distance(i->begin(),j)]=(*vMean)[distance(i->begin(),j)]+(*j)/(iNumOfRows-iBurnIn);

		}
	}


	for(vector<vector<double> >::iterator i=mData.begin(); i != mData.end(); ++i)
	{
		if(distance(mData.begin(),i)==iNextIndex-1)
		{
			for(vector<double>::iterator j= i->begin(); j != i -> end(); ++j)
			{
				(*mMatrix)[iCurrentMatrixIndex*iNumOfColumns+distance(i->begin(),j)]=*j;
			}
			iNextIndex=iNextIndex+iThin;
			iCurrentMatrixIndex=iCurrentMatrixIndex+1;
		}
	}


	iNumOfParam[0]=iNumOfMatrixRows;
	iNumOfVar[0]=iNumOfColumns;
	cout <<"The number of loaded rows is: "<<iNumOfRows<<endl;
	cout <<"The number of loaded columns is: "<<iNumOfColumns<<endl;
	cout <<"The number of imported rows is: "<<iNumOfMatrixRows<<endl;
	for(int j=0; j<iNumOfColumns;j++)
	{
		cout << "Mean "<<j<<" " << (*vMean)[j]<<endl;;
	}

}

void ImportData(string sFile,double** vData, double** vTimes, double** vPrice, double** vLogReturn,  int* iNumOfObs )
{
	vector<vector<double> > mData;

	fstream ifsDataFile;
	ifsDataFile.open(sFile.c_str(), ios::in);
	if(!ifsDataFile.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	string line;
	int iNumOfRows=0;
	while ( getline (ifsDataFile,line) )
	{
		vector<double>vRow;
		stringstream ssLine(line);
		string var;
		char* end;
		while(getline(ssLine, var,',') )
		{
			vRow.push_back(strtod( var.c_str() ,&end) );
		}
		mData.push_back(vRow);

		iNumOfRows=iNumOfRows+1;
		vRow.clear();
	}

	ifsDataFile.close();

	*vTimes = new double[iNumOfRows];
	*vData = new double[iNumOfRows];
	*vPrice = new double[iNumOfRows];
	*vLogReturn = new double[iNumOfRows];
	int k=0;
	for(vector<vector<double> >::iterator i=mData.begin(); i != mData.end(); ++i)
	{
//		cout<<(*i)[1]<<endl;
		(*vTimes)[k]= (*i)[0];
		(*vData)[k]= (*i)[1];
		(*vPrice)[k]= (*i)[2];
		(*vLogReturn)[k]= (*i)[3];
		k=k+1;
	}

	iNumOfObs[0]=iNumOfRows;
	cout <<"The number of observations is: "<<iNumOfRows<<endl;

}

void SetParametersSkellam(struct AllParam sParam, double * mParam, gsl_matrix * mWTilde , int iNumOfObs,  int iNumOfVar,  int iNumOfKnots, int iIter)
{
	sParam.dMu[0]=mParam[iNumOfVar*iIter];
	sParam.dPhi[0]=mParam[iNumOfVar*iIter+1];
	sParam.dSigma2[0]=mParam[iNumOfVar*iIter+2];
	sParam.dGamma[0]=mParam[iNumOfVar*iIter+3];

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+mParam[iNumOfVar*iIter+4+j]*gsl_matrix_get(mWTilde,i,j);
		}

	}


//	cout<<"XXXXXX Parameters XXXXX"<<endl;
//	cout << "Mu "<<sParam.dMu[0]<<endl;
//	cout << "Phi "<<sParam.dPhi[0]<<endl;
//	cout << "Sigma2 "<<sParam.dSigma2[0]<<endl;
//	cout << "Gamma "<<sParam.dGamma[0]<<endl;
//	cout << "Beta 1 "<< mParam[iNumOfVar*iIter+4]<<endl;
//	cout << "Beta 2 "<< mParam[iNumOfVar*iIter+5]<<endl;
//	cout << "Beta 3 "<< mParam[iNumOfVar*iIter+6]<<endl;
//	cout << "Beta 4 "<< mParam[iNumOfVar*iIter+7]<<endl;

}

void SetBICParametersSkellam(struct AllParam sParam, double * vMean, gsl_matrix * mWTilde , int iNumOfObs,  int iNumOfVar,  int iNumOfKnots)
{
	sParam.dMu[0]=vMean[0];
	sParam.dPhi[0]=vMean[1];
	sParam.dSigma2[0]=vMean[2];
	sParam.dGamma[0]=vMean[3];

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+vMean[4+j]*gsl_matrix_get(mWTilde,i,j);
		}

	}

}


void SetParametersDNB(struct AllParam sParam, double * mParam, gsl_matrix * mWTilde , int iNumOfObs,  int iNumOfVar,  int iNumOfKnots, int iIter)
{
	sParam.dMu[0]=mParam[iNumOfVar*iIter];
	sParam.dPhi[0]=mParam[iNumOfVar*iIter+1];
	sParam.dSigma2[0]=mParam[iNumOfVar*iIter+2];
	sParam.dGamma[0]=mParam[iNumOfVar*iIter+3];

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+mParam[iNumOfVar*iIter+4+j]*gsl_matrix_get(mWTilde,i,j);
		}

	}

	sParam.dNu[0]=mParam[iNumOfVar*iIter+iNumOfVar-1];
}

void SetBICParametersDNB(struct AllParam sParam, double * vMean, gsl_matrix * mWTilde , int iNumOfObs,  int iNumOfVar,  int iNumOfKnots)
{
	sParam.dMu[0]=vMean[0];
	sParam.dPhi[0]=vMean[1];
	sParam.dSigma2[0]=vMean[2];
	sParam.dGamma[0]=vMean[3];

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+vMean[4+j]*gsl_matrix_get(mWTilde,i,j);
		}

	}

	sParam.dNu[0]=vMean[iNumOfVar-1];
}


double DrawInital(double dPhi ,double dSigma2, const gsl_rng * gsl_random_num)
{
	return gsl_ran_gaussian(gsl_random_num,sqrt(dSigma2/(1-dPhi*dPhi)) );
}

double DrawTransition(double dPrevX, double dPhi ,double dSigma2, const gsl_rng * gsl_random_num)
{
	return dPhi * dPrevX+ gsl_ran_gaussian( gsl_random_num,sqrt(dSigma2)) ;
}

void SystemicResampling(double * vParticles, double * vCumNormalizedWeights,  int iNumOfParticles,
		const gsl_rng * gsl_random_num)
{

	double* vOldParticles=new double[iNumOfParticles];
	for(int i=0; i<iNumOfParticles; i++)
	{
		 vOldParticles[i]=vParticles[i];
	}

	double dU=gsl_rng_uniform(gsl_random_num);
	double dP;
	 int k=0;
	for(int i=0; i<iNumOfParticles; i++)
	{

		dP=dU/(double)iNumOfParticles+(double)i/(double)iNumOfParticles;
//		cout <<"dP " << dP << " vCumNormalizedWeights[ k]  "<<vCumNormalizedWeights[ k] <<endl;

		while(vCumNormalizedWeights[ k] <dP)
		{
			k=k+1;
		}

		vParticles[i]=vOldParticles[k];
	}

	delete [] vOldParticles;
}


double LogModifiedBesselFirstKindLargeXOut(int iN, double dX)
{
	int iNumOfK=200;
	double dValue=1;
	for(int i=iNumOfK; i>0;i--)
	{
		dValue=1-dValue*(4*pow((double) iN,2)-pow((double) 2*i-1,2))/(8*dX*i);
	}

	dValue=dX-0.5*log(2*PI*dX)+log(dValue);
	return dValue;
}

double ZeroSkellamPdfOut(const int iX, const double dGamma, const double dLambda){


	double dPdf;
	double d2Lambda;
	d2Lambda=2*dLambda;

	int iAbsX;
	iAbsX=abs(iX);
	gsl_sf_result result;


//	int status = gsl_sf_bessel_In_scaled_e(iAbsX,  d2Lambda ,&result);
	int status;
//	if(iAbsX>70)
//	{
//		status=0;
//		result.val=0;
//	}
//	else{
		gsl_set_error_handler_off();
		status =gsl_sf_bessel_In_scaled_e(iAbsX,  d2Lambda ,&result);
		gsl_set_error_handler (NULL);
//	}


	if (status) { /* an error occurred */
		/* status value specifies the type of error */
		double dApproxVal;
		if(iAbsX>70)
		{
			dApproxVal=0;
		}
		else
		{
			if(d2Lambda>2)
			{
				/*Large lambda approximation */
				dApproxVal=exp(-d2Lambda+LogModifiedBesselFirstKindLargeXOut(iAbsX, d2Lambda));
			}
			else
			{
				/*Small lambda approximation */
				dApproxVal=exp(-d2Lambda+iAbsX*log(dLambda)-lgamma(iAbsX+1));
			}
		}
		if(iX==0)
		{
			dPdf= dGamma + (1-dGamma)*dApproxVal;

		}
		else
		{
			dPdf=(1-dGamma)*dApproxVal;
		}
	}
	else
	{
		if(iX==0)
		{
			dPdf= dGamma + (1-dGamma)*result.val;

		}
		else
		{
			dPdf=(1-dGamma)*result.val;
		}
	}





	return dPdf;
}

void BootstrapFilterSkellamAdapt( struct AllParam sParam, double* vData,  int iNumOfObs,  int iNumOfParticles, double* vL)
{


//	double dMaxLogWeight=-1e20;
	double dSumWeights;
	double dESS=0;


	double* vParticles=new double[iNumOfParticles];
	double* vLogWeights=new double[iNumOfParticles];
	double* vWeights=new double[iNumOfParticles];
	double* vNormalizedWeights=new double[iNumOfParticles];
	double* vCumNormalizedWeights=new double[iNumOfParticles];

	double* vMeans=new double[iNumOfObs];

	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 2);



	for(int i=0; i<iNumOfObs; i++)
	{
		/* Initial step*/
		if(i==0)
		{
			dSumWeights=0;

			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Initital draw */
				vParticles[j]= DrawInital(sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "init particles " << vParticles[j] <<endl;

				/*Weights*/
				vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))/iNumOfParticles;

				/* Sum of weights */
				dSumWeights+=vWeights[j];
//				/* Log Weights*/
//				vLogWeights[j]=log(ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j])));
//
//				/* max of log weights */
//				if(dMaxLogWeight<vLogWeights[j])
//				{
//					dMaxLogWeight=vLogWeights[j];
//				}
			}



			/* Calculate statistics */
			dESS=0;
			vL[i]=0;
			vMeans[i]=0;
			/* Normalized weights */
			for(int j=0; j<iNumOfParticles; j++)
			{
				/*Calculate likelihood*/
				vL[i]=vL[i]+vWeights[j];

				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
//				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights " << vNormalizedWeights[j]<<endl;
				if(j==0)
				{
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights " << vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;

//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
		}
		else
		{
//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
			/*Resample */
			if(dESS<0.5*iNumOfParticles)
			{

//				cout << "XXX Resampling  XXX"<<endl;
				//cout << "XXX ESS "<<dESS<<" XXX"<<endl;

//				for(int k=0; k<iNumOfParticles; k++)
//				{
//					cout<<"Old Partciles " << vParticles[k] <<" Weights " << vNormalizedWeights[k] <<endl;
//				}
				SystemicResampling(vParticles, vCumNormalizedWeights,iNumOfParticles,gsl_random_num);
//				for(int k=0; k<iNumOfParticles; k++)
//				{
//					cout<<"New Partciles " << vParticles[k] <<endl;
//				}
			}

			dSumWeights=0;
			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Propagate */
				vParticles[j]=DrawTransition(vParticles[j],sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "particles " << vParticles[j] <<endl;
				/*Weights*/
				if(dESS<0.5*iNumOfParticles)
				{
					vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))/iNumOfParticles;
				}
				else
				{
					vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))*vNormalizedWeights[j];
				}
//				cout << "weights " << vWeights[j]<<endl;
				/* Sum of weights */
				dSumWeights+=vWeights[j];
			}






			dESS=0;
			vMeans[i]=0;
			vL[i]=0;
			for(int j=0; j<iNumOfParticles; j++)
			{

				/*Calculate likelihood*/
				vL[i]=vL[i]+vWeights[j];

				/* Normalized weights */
				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				/* Calculate statistics */
//				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights at " << j <<" is "<< vNormalizedWeights[j]<<endl;
				if(j==0)
				{
//					cout<< "first "<<vNormalizedWeights[0] <<endl;
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
//					cout << "prev cum normalized weights at " << j-1 <<" is "<<  vCumNormalizedWeights[j-1]<<endl;
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights at " << j <<" is "<<  vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;
		}


	}



	delete [] vParticles;
	delete [] vMeans;
	delete [] vLogWeights;
	delete [] vNormalizedWeights;
	delete [] vWeights;

}

void BootstrapFilterSkellam( struct AllParam sParam, double* vData,  int iNumOfObs,  int iNumOfParticles, double* vL)
{


//	double dMaxLogWeight=-1e20;
	double dSumWeights;
	double dESS;


	double* vParticles=new double[iNumOfParticles];
	double* vLogWeights=new double[iNumOfParticles];
	double* vWeights=new double[iNumOfParticles];
	double* vNormalizedWeights=new double[iNumOfParticles];
	double* vCumNormalizedWeights=new double[iNumOfParticles];

	double* vMeans=new double[iNumOfObs];

	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 2);



	for(int i=0; i<iNumOfObs; i++)
	{
		/* Initial step*/
		if(i==0)
		{
			dSumWeights=0;

			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Initital draw */
				vParticles[j]= DrawInital(sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "init particles " << vParticles[j] <<endl;

				/*Weights*/
				vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))/iNumOfParticles;

				/* Sum of weights */
				dSumWeights+=vWeights[j];
//				/* Log Weights*/
//				vLogWeights[j]=log(ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j])));
//
//				/* max of log weights */
//				if(dMaxLogWeight<vLogWeights[j])
//				{
//					dMaxLogWeight=vLogWeights[j];
//				}
			}



			/* Calculate statistics */
			dESS=0;
			vL[i]=0;
			vMeans[i]=0;
			/* Normalized weights */
			for(int j=0; j<iNumOfParticles; j++)
			{
				/*Calculate likelihood*/
				vL[i]=vL[i]+vWeights[j];

				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights " << vNormalizedWeights[j]<<endl;
				if(j==0)
				{
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights " << vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;

//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
		}
		else
		{
//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
			/*Resample */

			SystemicResampling(vParticles, vCumNormalizedWeights,iNumOfParticles,gsl_random_num);
//

			dSumWeights=0;
			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Propagate */
				vParticles[j]=DrawTransition(vParticles[j],sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "particles " << vParticles[j] <<endl;
				/*Weights*/

				vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))/iNumOfParticles;

//				cout << "weights " << vWeights[j]<<endl;
				/* Sum of weights */
				dSumWeights+=vWeights[j];
			}






			dESS=0;
			vMeans[i]=0;
			vL[i]=0;
			for(int j=0; j<iNumOfParticles; j++)
			{

				/*Calculate likelihood*/
				vL[i]=vL[i]+vWeights[j];

				/* Normalized weights */
				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				/* Calculate statistics */
				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights at " << j <<" is "<< vNormalizedWeights[j]<<endl;
				if(j==0)
				{
//					cout<< "first "<<vNormalizedWeights[0] <<endl;
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
//					cout << "prev cum normalized weights at " << j-1 <<" is "<<  vCumNormalizedWeights[j-1]<<endl;
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights at " << j <<" is "<<  vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;
		}


	}

	string sPFFile="vPFMean.csv";
	WriteOutDoubleArrayOut(vMeans, iNumOfObs, 1, sPFFile);

	delete [] vParticles;
	delete [] vMeans;
	delete [] vLogWeights;
	delete [] vNormalizedWeights;
	delete [] vWeights;

}

double ZeroDNBPdfOut(const int iX, const double dGamma, const double dLambda, const double dNu)
{

		double dPdf;
		double dLambdaRatio=dLambda/(dLambda+dNu);
		double dNuRatio=dNu/(dLambda+dNu);

		int iAbsX;
		iAbsX=abs(iX);

		gsl_sf_result result;
		gsl_set_error_handler_off();
		int status = gsl_sf_hyperg_2F1_e(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) ,&result);
		gsl_set_error_handler (NULL);

		if (status) { /* an error occurred */
			/* status value specifies the type of error */

			if(iX==0)
			{
				dPdf= dGamma + (1-dGamma)*pow(dNuRatio,2*dNu)*HyperGeoOut(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) );
			}
			else
			{
				dPdf=(1-dGamma)*pow(dNuRatio,2*dNu)*pow(dLambdaRatio,iAbsX)*exp(gsl_sf_lngamma(dNu+iAbsX)-gsl_sf_lngamma(iAbsX+1)-gsl_sf_lngamma(dNu))*HyperGeoOut(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) );
			}
		}
		else{

			if(iX==0)
			{
				dPdf= dGamma + (1-dGamma)*pow(dNuRatio,2*dNu)*result.val;
			}
			else
			{
				dPdf=(1-dGamma)*pow(dNuRatio,2*dNu)*pow(dLambdaRatio,iAbsX)*exp(gsl_sf_lngamma(dNu+iAbsX)-gsl_sf_lngamma(iAbsX+1)-gsl_sf_lngamma(dNu))*result.val;
			}
		}



////		cout<<"iX "<<iX<<endl;
////		cout<<"dGamma "<<dGamma<<endl;
////		cout<<"dNuRatio "<<dNuRatio<<endl;
////		cout<<"dLambdaRatio "<<dLambdaRatio<<endl;
////		cout<<"dNuRatio "<<dNuRatio<<endl;
////
////		cout<<"XXXXXX Check HyperGeo XXXXXXXXXXXX"<<endl;
////		double dGSL=gsl_sf_hyperg_2F1(dNu+abs(iX),dNu, abs(iX)+1 , pow(dLambdaRatio,2) );
////		double dTaylor=HyperGeoOut(dNu+abs(iX),dNu, abs(iX)+1 , pow(dLambdaRatio,2) );
////		cout<<"GSL -> "<<dGSL <<endl;
////		cout<<"Taylor ->  "<<dTaylor <<endl;
////		cout<<"Diff -> "<< dTaylor-dGSL<<endl;
//		if(iX==0)
//		{
//			dPdf= dGamma + (1-dGamma)*pow(dNuRatio,2*dNu)*gsl_sf_hyperg_2F1(dNu,dNu, 1 , pow(dLambdaRatio,2));
//		}
//		else
//		{
//			int iAbsX;
//			iAbsX=abs(iX);
////			dPdf=(1-dGamma)*pow(dNuRatio,2*dNu)*pow(dLambdaRatio,iAbsX)*(gsl_sf_gamma(dNu+iAbsX)/(gsl_sf_gamma(iAbsX+1)*gsl_sf_gamma(dNu)))*gsl_sf_hyperg_2F1(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) );
//
//
//
//			gsl_sf_result result;
//			int status = gsl_sf_hyperg_2F1_e(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) ,&result);
//			if (status) { /* an error occurred */
//
//				/* status value specifies the type of error */
//				dPdf=(1-dGamma)*pow(dNuRatio,2*dNu)*pow(dLambdaRatio,iAbsX)*exp(gsl_sf_lngamma(dNu+iAbsX)-gsl_sf_lngamma(iAbsX+1)-gsl_sf_lngamma(dNu))*HyperGeoOut(dNu+iAbsX,dNu, iAbsX+1 , pow(dLambdaRatio,2) );
//			}
//			else{
//				dPdf=(1-dGamma)*pow(dNuRatio,2*dNu)*pow(dLambdaRatio,iAbsX)*exp(gsl_sf_lngamma(dNu+iAbsX)-gsl_sf_lngamma(iAbsX+1)-gsl_sf_lngamma(dNu))*result.val;
//			}
//
//		}


		return dPdf;
}

void BootstrapFilterDNB( struct AllParam sParam, double* vData,  int iNumOfObs,  int iNumOfParticles, double* vL)
{


//	double dMaxLogWeight=-1e20;
	double dSumWeights;
	double dESS;


	double* vParticles=new double[iNumOfParticles];
	double* vLogWeights=new double[iNumOfParticles];
	double* vWeights=new double[iNumOfParticles];
	double* vNormalizedWeights=new double[iNumOfParticles];
	double* vCumNormalizedWeights=new double[iNumOfParticles];

	double* vMeans=new double[iNumOfObs];

	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 2);



	for(int i=0; i<iNumOfObs; i++)
	{
		/* Initial step*/
		if(i==0)
		{
			dSumWeights=0;

			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Initital draw */
				vParticles[j]= DrawInital(sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "init particles " << vParticles[j] <<endl;

				/*Weights*/
				vWeights[j]=ZeroDNBPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]), sParam.dNu[0])/iNumOfParticles;

				/* Sum of weights */
				dSumWeights+=vWeights[j];
//				/* Log Weights*/
//				vLogWeights[j]=log(ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j])));
//
//				/* max of log weights */
//				if(dMaxLogWeight<vLogWeights[j])
//				{
//					dMaxLogWeight=vLogWeights[j];
//				}
			}



			/* Calculate statistics */
			dESS=0;
			vL[i]=0;
			vMeans[i]=0;
			/* Normalized weights */
			for(int j=0; j<iNumOfParticles; j++)
			{
				/*Calculate likelihood*/
				vL[i]=vL[i]+vWeights[j];

				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights " << vNormalizedWeights[j]<<endl;
				if(j==0)
				{
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights " << vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;

//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
		}
		else
		{
//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
			/*Resample */

			SystemicResampling(vParticles, vCumNormalizedWeights,iNumOfParticles,gsl_random_num);
//

			dSumWeights=0;
			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Propagate */
				vParticles[j]=DrawTransition(vParticles[j],sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "particles " << vParticles[j] <<endl;
				/*Weights*/

				vWeights[j]=ZeroDNBPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]), sParam.dNu[0])/iNumOfParticles;

//				cout << "weights " << vWeights[j]<<endl;
				/* Sum of weights */
				dSumWeights+=vWeights[j];
			}






			dESS=0;
			vMeans[i]=0;
			vL[i]=0;
			for(int j=0; j<iNumOfParticles; j++)
			{

				/*Calculate likelihood*/
				vL[i]=vL[i]+vWeights[j];

				/* Normalized weights */
				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				/* Calculate statistics */
				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights at " << j <<" is "<< vNormalizedWeights[j]<<endl;
				if(j==0)
				{
//					cout<< "first "<<vNormalizedWeights[0] <<endl;
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
//					cout << "prev cum normalized weights at " << j-1 <<" is "<<  vCumNormalizedWeights[j-1]<<endl;
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights at " << j <<" is "<<  vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;
		}


	}

	string sPFFile="vPFMean.csv";
	WriteOutDoubleArrayOut(vMeans, iNumOfObs, 1, sPFFile);

	delete [] vParticles;
	delete [] vMeans;
	delete [] vLogWeights;
	delete [] vNormalizedWeights;
	delete [] vWeights;

}

void OneStepAheadSkellam( struct AllParam sParam, double* vData, double dLow, double dHigh,  double* vPrice, int iNumOfObs,  int iNumOfParticles,
		double* vPredMean, double* mPredDist)
{

	double dRange=dHigh-dLow;

//	double dMaxLogWeight=-1e20;
	double dSumWeights;
	double dESS;

	for(int i=0; i<iNumOfObs; i++)
	{
		for(int k=0; k<int(dRange)+1 ; k++)
		{
			mPredDist[i*(int(dRange)+1)+k]=0;
		}
	}


	for(int k=0; k<int(dRange)+1 ; k++)
	{
		mPredDist[k]=dLow+k;
	}

	double* vParticles=new double[iNumOfParticles];
//	double* vForecast=new double[iNumOfParticles];
	double* vLogWeights=new double[iNumOfParticles];
	double* vWeights=new double[iNumOfParticles];
	double* vNormalizedWeights=new double[iNumOfParticles];
	double* vCumNormalizedWeights=new double[iNumOfParticles];

	double* vMeansVol=new double[iNumOfObs];

	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 2);



	for(int i=0; i<iNumOfObs; i++)
	{
		cout<<"Time "<<i<<endl;
//		cout<<"data "<<vData[i]<<endl;
		/* Initial step*/
		if(i==0)
		{
			dSumWeights=0;

			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Initital draw */
				vParticles[j]= DrawInital(sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "init particles " << vParticles[j] <<endl;

				/*Weights*/
				vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))/iNumOfParticles;

				/* Sum of weights */
				dSumWeights+=vWeights[j];
//				/* Log Weights*/
//				vLogWeights[j]=log(ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j])));
//
//				/* max of log weights */
//				if(dMaxLogWeight<vLogWeights[j])
//				{
//					dMaxLogWeight=vLogWeights[j];
//				}
			}



			/* Calculate statistics */
			dESS=0;
			vPredMean[i]=0;
			vMeansVol[i]=0;
			/* Normalized weights */
			for(int j=0; j<iNumOfParticles; j++)
			{
				/*Calculate likelihood*/


				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				vMeansVol[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights " << vNormalizedWeights[j]<<endl;
				if(j==0)
				{
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights " << vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;

//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
		}
		else
		{
//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
			/*Resample */

			SystemicResampling(vParticles, vCumNormalizedWeights,iNumOfParticles,gsl_random_num);
//

			dSumWeights=0;
			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Propagate */
				vParticles[j]=DrawTransition(vParticles[j],sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);

//				cout << "particles " << vParticles[j] <<endl;
				/*Weights*/

				vWeights[j]=ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]))/iNumOfParticles;

//				cout << "weights " << vWeights[j]<<endl;
				/* Sum of weights */
				dSumWeights+=vWeights[j];
			}





			double dProb=0;
			double dTempProb;
			double dThreshold=0.000000001;
			for(int j=0; j<iNumOfParticles; j++)
			{
				for(int k=0; k<(int(dRange)+1); k++)
				{
					double dCurrentValue=dLow+k;


					if(k==0)
					{

						dProb=0;
						dTempProb=1;
						int m=0;
						while((dTempProb>dThreshold)&(m<100))
						{
							dTempProb=ZeroSkellamPdfOut(dCurrentValue-m,  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i-1]+vParticles[j]));
							dProb+=dTempProb;
							m+=1;
						}
					}
					else if(k==int(dRange))
					{
						dProb=0;
						dTempProb=1;
						int m=0;
						while((dTempProb>dThreshold)&(m<100))
						{
							dTempProb=ZeroSkellamPdfOut(dCurrentValue+m,  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i-1]+vParticles[j]));
							dProb+=dTempProb;
							m+=1;
						}

					}
					else
					{
						dProb=ZeroSkellamPdfOut(dCurrentValue,  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i-1]+vParticles[j]));
					}


					mPredDist[i*(int(dRange)+1)+k]+=dProb* vNormalizedWeights[j];
				}
			}




			dESS=0;
			vMeansVol[i]=0;
			vPredMean[i]=0;
			for(int j=0; j<iNumOfParticles; j++)
			{


				/* Normalized weights */
				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				/* Calculate statistics */
				vMeansVol[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights at " << j <<" is "<< vNormalizedWeights[j]<<endl;
				if(j==0)
				{
//					cout<< "first "<<vNormalizedWeights[0] <<endl;
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
//					cout << "prev cum normalized weights at " << j-1 <<" is "<<  vCumNormalizedWeights[j-1]<<endl;
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights at " << j <<" is "<<  vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;
		}


	}

	string sVolFile="vVolMean.csv";
	WriteOutDoubleArrayOut(vMeansVol, iNumOfObs, 1, sVolFile);

//	string sPredMeanFile="vPredMean.csv";
//	WriteOutDoubleArrayOut(vPredMean, iNumOfObs, 1, sPredMeanFile);

	delete [] vParticles;
	delete [] vMeansVol;
	delete [] vLogWeights;
	delete [] vNormalizedWeights;
	delete [] vWeights;

}

void OneStepAheadDNB( struct AllParam sParam, double* vData, double dLow, double dHigh,  double* vPrice, int iNumOfObs,  int iNumOfParticles,
		double* vPredMean, double* mPredDist)
{

	double dRange=dHigh-dLow;

//	double dMaxLogWeight=-1e20;
	double dSumWeights;
	double dESS;

	for(int i=0; i<iNumOfObs; i++)
	{
		for(int k=0; k<int(dRange)+1 ; k++)
		{
			mPredDist[i*(int(dRange)+1)+k]=0;
		}
	}

	for(int k=0; k<int(dRange)+1 ; k++)
	{
		mPredDist[k]=dLow+k;
	}



	double* vParticles=new double[iNumOfParticles];
//	double* vForecast=new double[iNumOfParticles];
	double* vLogWeights=new double[iNumOfParticles];
	double* vWeights=new double[iNumOfParticles];
	double* vNormalizedWeights=new double[iNumOfParticles];
	double* vCumNormalizedWeights=new double[iNumOfParticles];

	double* vMeansVol=new double[iNumOfObs];

	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 2);



	for(int i=0; i<iNumOfObs; i++)
	{
		cout<<"Time "<<i<<endl;
//		cout<<"Low "<<dLow<<endl;
//		cout<<"High "<<dHigh<<endl;
		/* Initial step*/
		if(i==0)
		{

			dSumWeights=0;

			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Initital draw */
				vParticles[j]= DrawInital(sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);
//				cout << "init particles " << vParticles[j] <<endl;

				/*Weights*/
				vWeights[j]=ZeroDNBPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]), sParam.dNu[0])/iNumOfParticles;

				/* Sum of weights */
				dSumWeights+=vWeights[j];
//				/* Log Weights*/
//				vLogWeights[j]=log(ZeroSkellamPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j])));
//
//				/* max of log weights */
//				if(dMaxLogWeight<vLogWeights[j])
//				{
//					dMaxLogWeight=vLogWeights[j];
//				}
			}



			/* Calculate statistics */
			dESS=0;
			vPredMean[i]=0;
			vMeansVol[i]=0;
			/* Normalized weights */
			for(int j=0; j<iNumOfParticles; j++)
			{
				/*Calculate likelihood*/


				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				vMeansVol[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights " << vNormalizedWeights[j]<<endl;
				if(j==0)
				{
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights " << vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;

//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
		}
		else
		{
//			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
			/*Resample */

			SystemicResampling(vParticles, vCumNormalizedWeights,iNumOfParticles,gsl_random_num);
//

			dSumWeights=0;
			for(int j=0; j<iNumOfParticles; j++)
			{
				/* Propagate */
				vParticles[j]=DrawTransition(vParticles[j],sParam.dPhi[0] ,sParam.dSigma2[0],  gsl_random_num);

//				cout << "particles " << vParticles[j] <<endl;
				/*Weights*/

				vWeights[j]=ZeroDNBPdfOut(vData[i],  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i]+vParticles[j]), sParam.dNu[0])/iNumOfParticles;

//				cout << "weights " << vWeights[j]<<endl;
				/* Sum of weights */
				dSumWeights+=vWeights[j];
			}





			double dProb=0;
			double dTempProb;
			double dThreshold=0.000000001;
			for(int k=0; k<(int(dRange)+1); k++)
			{
				for(int j=0; j<iNumOfParticles; j++)
				{

					double dCurrentValue=dLow+k;
//					cout<<"Time "<<i<<endl;
//					cout<<"Data "<<vData[i]<<endl;
//					cout<<"dCurrent "<<dCurrentValue<<endl;

					if(k==0)
					{

						dProb=0;
						dTempProb=1;
						int m=0;
						while((dTempProb>dThreshold)&(m<100))
						{
							dTempProb=ZeroDNBPdfOut(dCurrentValue-m,  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i-1]+vParticles[j]), sParam.dNu[0]);
							dProb+=dTempProb;
							m+=1;
						}
					}
					else if(k==int(dRange))
					{
						dProb=0;
						dTempProb=1;
						int m=0;
						while((dTempProb>dThreshold)&(m<100))
						{
							dTempProb=ZeroDNBPdfOut(dCurrentValue+m,  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i-1]+vParticles[j]), sParam.dNu[0]);
							dProb+=dTempProb;
							m+=1;
						}

					}
					else
					{
						dProb=ZeroDNBPdfOut(dCurrentValue,  sParam.dGamma[0], exp(sParam.dMu[0]+sParam.vS[i-1]+vParticles[j]), sParam.dNu[0]);
					}

//					cout<<"k is "<< k<<endl;
//					cout<<"dProb is "<<dProb<<endl;
//					cout<<"Weight is "<<vNormalizedWeights[j]<<endl;
//					cout<<"Dist cont  is "<<dProb* vNormalizedWeights[j]<<endl;
//					cout<<"Dist is "<<mPredDist[i*(int(dRange)+1)+k]<<endl;
					mPredDist[i*(int(dRange)+1)+k]+=dProb* vNormalizedWeights[j];
				}
			}

//			for(int k=0; k<(int(dRange)+1); k++)
//			{
//				cout<<"XXXX Dist is "<<mPredDist[i*(int(dRange)+1)+k]<<endl;
//			}


			dESS=0;
			vMeansVol[i]=0;
			vPredMean[i]=0;
			for(int j=0; j<iNumOfParticles; j++)
			{


				/* Normalized weights */
				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
				/* Calculate statistics */
				vMeansVol[i]+=vNormalizedWeights[j]*vParticles[j];
//				cout << "normalized weights at " << j <<" is "<< vNormalizedWeights[j]<<endl;
				if(j==0)
				{
//					cout<< "first "<<vNormalizedWeights[0] <<endl;
					vCumNormalizedWeights[j]=vNormalizedWeights[j];
				}
				else
				{
//					cout << "prev cum normalized weights at " << j-1 <<" is "<<  vCumNormalizedWeights[j-1]<<endl;
					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
				}
//				cout << "cum normalized weights at " << j <<" is "<<  vCumNormalizedWeights[j]<<endl;
				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
			}
			dESS=1/dESS;
		}


	}

	string sVolFile="vVolMean.csv";
	WriteOutDoubleArrayOut(vMeansVol, iNumOfObs, 1, sVolFile);

//	string sPredMeanFile="vPredMean.csv";
//	WriteOutDoubleArrayOut(vPredMean, iNumOfObs, 1, sPredMeanFile);

	delete [] vParticles;
	delete [] vMeansVol;
	delete [] vLogWeights;
	delete [] vNormalizedWeights;
	delete [] vWeights;

}

#undef PI
#endif /* BASIC_FUNCTION2_H_ */
