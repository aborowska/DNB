/*
 * OutOfSample.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: istvan
 */

#include "OutOfSample.h"
#include "basic_functions2.h"
#include <time.h>

OutOfSample::OutOfSample(string sFileParam,unsigned int iBurnIn, unsigned int iThin, string sType, string sFileDataIn,string sFileDataOut,string sPath)
{
	/* Import parameter draws and thin */

	CreatPrefixOut(sFileDataIn,sType, &sPrefix );
	cout<< "prefix done" <<endl;
	ImportMatrix(sPrefix+sFileParam,  iBurnIn, iThin, &mParam, &vMean, &iNumOfParam,&iNumOfVar);
	cout<< "Check" <<endl;
	cout<< sPrefix+sFileParam<<endl;
	/* Import data */
	// ImportData(sFileDataIn,&vDataIn, &vTimesIn,&vPriceIn, &vLogReturnIn, &iNumOfObsIn);
	ImportData2(sPath,sFileDataIn,&vDataIn, &vTimesIn,&vPriceIn, &vLogReturnIn, &iNumOfObsIn);
	/* Import data */
	// ImportData(sFileDataOut,&vDataOut, &vTimesOut,&vPriceOut, &vLogReturnOut,&iNumOfObsOut);
	ImportData2(sPath,sFileDataOut,&vDataOut, &vTimesOut,&vPriceOut, &vLogReturnOut,&iNumOfObsOut);
// private:

	// /* Number of observations */
	// int iNumOfObsIn;
	// int iNumOfObsOut;

	// /* Number of parameters */
	// int iNumOfParam;

	// /* Number of variables */
	// int iNumOfVar;
	// /* Data */
	// double* vDataIn;
	// double* vTimesIn;
	// double* vPriceIn;
	// double* vLogReturnIn;
	// double* vDataOut;
	// double* vTimesOut;
	// double* vPriceOut;
	// double* vLogReturnOut;
	// double* mParam;
	// double* vMean;
	// string sPrefix;	

}




void OutOfSample::GoodnessOfFitOrdNormal(int iNumOfParticles)
{

	/* Calculate W matrix for the seasonality*/
	int iNumOfKnots=3;
	double * vKnots= new double[iNumOfKnots];
	double * vTempL=new double[iNumOfObsOut];
	double * vLogL=new double[iNumOfObsOut];
//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=5400; /* 11 */
//	vKnots[2]=10800; /* 12.30 */
//	vKnots[3]=18000;/* 14.30 */
//	vKnots[iNumOfKnots-1]=23100; /* 15.55 */

	vKnots[0]=0;
	vKnots[1]=10800;
	vKnots[iNumOfKnots-1]=23400;

	gsl_matrix * mWTildeIn = gsl_matrix_alloc (iNumOfObsIn, iNumOfKnots-1);
	gsl_matrix * mWTildeOut = gsl_matrix_alloc (iNumOfObsOut, iNumOfKnots-1);
	CalculateSplineMatricesOut(iNumOfObsIn,  vTimesIn, iNumOfKnots,  vKnots,  mWTildeIn );
	CalculateSplineMatricesOut(iNumOfObsOut,  vTimesOut, iNumOfKnots,  vKnots,  mWTildeOut );

	/* Parameters */
	AllParam sParam;
	sParam.dMu=new double;
	sParam.dPhi=new double;
	sParam.dSigma2=new double;
	sParam.dGamma=new double;
	sParam.vS=new double [iNumOfObsOut];

	for(int i=0;i<iNumOfObsOut; i++)
	{
		vLogL[i]=0;
	}
	/* For loop over the thinned parameters */
	for(int i=0;i<iNumOfParam; i++)
	{
		cout << "Param iteration: "<< i<<endl;
		/* Set parameters to param struct */
		SetParametersNormal( sParam, mParam, mWTildeOut ,iNumOfObsOut, iNumOfVar,  iNumOfKnots, i);

		/* Run particle filter */
		time_t start,end;
		time (&start);
		//BootstrapFilterSkellamAdapt( sParam,  vData, iNumOfObs, iNumOfParticles,  vTempL);
		BootstrapFilterOrderedNormal( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
		time (&end);
		double dif = difftime (end,start);
		printf ("Elasped time is %.2lf seconds. \n", dif );
		/* Average results */
		for(int j=0;j<iNumOfObsOut; j++)
		{
			vLogL[j]=vLogL[j]+vTempL[j]/iNumOfParam;
		}

	}


	for(int j=0;j<iNumOfObsOut; j++)
	{
		vLogL[j]=log(vLogL[j]);
	}


	string sLLFile="vLogLikeNormal.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sLLFile);
	cout << "Predictive LL done!"<<endl;

	/* BIC Out of sample */
	SetBICParametersNormal(sParam,  vMean,  mWTildeOut , iNumOfObsOut,   iNumOfVar, iNumOfKnots);
	BootstrapFilterOrderedNormal( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
	for(int j=0;j<iNumOfObsOut; j++)
	{
			vLogL[j]=log(vTempL[j]);
	}

	string sOutBICFile="vOutNormalBIC.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sOutBICFile);
	cout << "BIC Out done!"<<endl;

	/* BIC In sample */

	double * vTempLIn=new double[iNumOfObsIn];
	double * vLogLIn=new double[iNumOfObsIn];

	delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsIn];

	SetBICParametersNormal(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
	string sInSeasonal="vS_Normal_Check.csv";
	WriteOutDoubleArrayOut(sParam.vS, iNumOfObsIn, 1, sPrefix+sInSeasonal);

	BootstrapFilterOrderedNormal( sParam,  vDataIn, iNumOfObsIn, iNumOfParticles,  vTempLIn);
	for(int j=0;j<iNumOfObsIn; j++)
	{
		vLogLIn[j]=log(vTempLIn[j]);
	}

	string sInBICFile="vInNormalBIC.csv";
	WriteOutDoubleArrayOut(vLogLIn, iNumOfObsIn, 1, sPrefix+sInBICFile);

	cout << "BIC In done!"<<endl;



	delete [] vTempL;
	delete [] vLogL;
	delete [] vTempLIn;
	delete [] vLogLIn;
	delete [] vKnots;

	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;

	delete [] sParam.vS;
}



void OutOfSample::GoodnessOfFitOrdT(int iNumOfParticles)
{

	/* Calculate W matrix for the seasonality*/
	int iNumOfKnots=3;
	double * vKnots= new double[iNumOfKnots];
	double * vTempL=new double[iNumOfObsOut];
	double * vLogL=new double[iNumOfObsOut];
//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=5400; /* 11 */
//	vKnots[2]=10800; /* 12.30 */
//	vKnots[3]=18000;/* 14.30 */
//	vKnots[iNumOfKnots-1]=23100; /* 15.55 */

	vKnots[0]=0;
	vKnots[1]=10800;
	vKnots[iNumOfKnots-1]=23400;

	gsl_matrix * mWTildeIn = gsl_matrix_alloc (iNumOfObsIn, iNumOfKnots-1);
	gsl_matrix * mWTildeOut = gsl_matrix_alloc (iNumOfObsOut, iNumOfKnots-1);
	CalculateSplineMatricesOut(iNumOfObsIn,  vTimesIn, iNumOfKnots,  vKnots,  mWTildeIn );
	CalculateSplineMatricesOut(iNumOfObsOut,  vTimesOut, iNumOfKnots,  vKnots,  mWTildeOut );

	/* Parameters */
	AllParam sParam;
	sParam.dMu=new double;
	sParam.dPhi=new double;
	sParam.dSigma2=new double;
	sParam.dGamma=new double;
	sParam.dNu=new double;
	sParam.vS=new double [iNumOfObsOut];

	for(int i=0;i<iNumOfObsOut; i++)
	{
		vLogL[i]=0;
	}
	/* For loop over the thinned parameters */
	for(int i=0;i<iNumOfParam; i++)
	{
		cout << "Param iteration: "<< i<<endl;
		/* Set parameters to param struct */
		SetParametersT( sParam, mParam, mWTildeOut ,iNumOfObsOut, iNumOfVar,  iNumOfKnots, i);

		/* Run particle filter */
		time_t start,end;
		time (&start);
		//BootstrapFilterSkellamAdapt( sParam,  vData, iNumOfObs, iNumOfParticles,  vTempL);
		BootstrapFilterOrderedT( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
		time (&end);
		double dif = difftime (end,start);
		printf ("Elasped time is %.2lf seconds. \n", dif );
		/* Average results */
		for(int j=0;j<iNumOfObsOut; j++)
		{
			vLogL[j]=vLogL[j]+vTempL[j]/iNumOfParam;
		}

	}


	for(int j=0;j<iNumOfObsOut; j++)
	{
		vLogL[j]=log(vLogL[j]);
	}


	string sLLFile="vLogLikeT.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sLLFile);
	cout << "Predictive LL done!"<<endl;

	/* BIC Out of sample */
	SetBICParametersT(sParam,  vMean,  mWTildeOut , iNumOfObsOut,   iNumOfVar, iNumOfKnots);
	BootstrapFilterOrderedT( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
	for(int j=0;j<iNumOfObsOut; j++)
	{
			vLogL[j]=log(vTempL[j]);
	}

	string sOutBICFile="vOutTBIC.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sOutBICFile);
	cout << "BIC Out done!"<<endl;

	/* BIC In sample */

	double * vTempLIn=new double[iNumOfObsIn];
	double * vLogLIn=new double[iNumOfObsIn];

	//delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsIn];

	SetBICParametersT(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
	BootstrapFilterOrderedT( sParam,  vDataIn, iNumOfObsIn, iNumOfParticles,  vTempLIn);
	for(int j=0;j<iNumOfObsIn; j++)
	{
		vLogLIn[j]=log(vTempLIn[j]);
	}

	string sInBICFile="vTBIC.csv";
	WriteOutDoubleArrayOut(vLogLIn, iNumOfObsIn, 1, sPrefix+sInBICFile);

	cout << "BIC In done!"<<endl;



	delete [] vTempL;
	delete [] vLogL;
	delete [] vTempLIn;
	delete [] vLogLIn;
	delete [] vKnots;

	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;
	delete sParam.dNu;

	delete [] sParam.vS;
}

void OutOfSample::ForecastOrdNormal(int iNumOfParticles)
{

	/* Calculate W matrix for the seasonality*/
	int iNumOfKnots=3;
	double * vKnots= new double[iNumOfKnots];


//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=5400; /* 11 */
//	vKnots[2]=10800; /* 12.30 */
//	vKnots[3]=18000;/* 14.30 */
//	vKnots[iNumOfKnots-1]=23100; /* 15.55 */

	vKnots[0]=0;
	vKnots[1]=10800;
	vKnots[iNumOfKnots-1]=23400;

	gsl_matrix * mWTildeIn = gsl_matrix_alloc (iNumOfObsIn, iNumOfKnots-1);
	gsl_matrix * mWTildeOut = gsl_matrix_alloc (iNumOfObsOut, iNumOfKnots-1);
	CalculateSplineMatricesOut(iNumOfObsIn,  vTimesIn, iNumOfKnots,  vKnots,  mWTildeIn );
	CalculateSplineMatricesOut(iNumOfObsOut,  vTimesOut, iNumOfKnots,  vKnots,  mWTildeOut );

	/* Parameters */
	AllParam sParam;
	sParam.dMu=new double;
	sParam.dPhi=new double;
	sParam.dSigma2=new double;
	sParam.dGamma=new double;
	sParam.vS=new double [iNumOfObsIn];

	SetBICParametersNormal(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
	double dLow=0;
	double dHigh=0;

	for(int i=0; i<iNumOfObsIn;i++)
	{
		if(vDataIn[i]<dLow)
		{
			dLow=vDataIn[i];
		}
		if(vDataIn[i]>dHigh)
		{
			dHigh=vDataIn[i];
		}
	}

	cout<<"dLow "<<dLow<<endl;
	cout<<"dHigh "<<dHigh<<endl;

	if(dHigh<dLow)
	{
		cout<<"!!!!!!!! Error dHigh < dLow !!!!!!!!!"<<endl;
	}
	double * vPredMeanIn=new double[iNumOfObsIn];
	double*  mPredDistIn=new double[iNumOfObsIn*(int(dHigh-dLow)+1)];
	BootstrapFilterOrderedNormalForecast(sParam, vDataIn,  dLow, dHigh,  iNumOfObsIn,  iNumOfParticles,vPredMeanIn,  mPredDistIn);

	string sPredDistFile="InPredDist.csv";
	WriteOutDoubleArrayOut(mPredDistIn, iNumOfObsIn, ((dHigh-dLow)+1), sPrefix+sPredDistFile);
	string sPredMeanFile="InPredMean.csv";
	WriteOutDoubleArrayOut(vPredMeanIn, iNumOfObsIn, 1, sPrefix+sPredMeanFile);


	cout << "Forecast In done!"<<endl;
	dLow=0;
	dHigh=0;
	for(int i=0; i<iNumOfObsOut;i++)
	{
		if(vDataOut[i]<dLow)
		{
			dLow=vDataOut[i];
		}
		if(vDataOut[i]>dHigh)
		{
			dHigh=vDataOut[i];
		}
	}

	cout<<"dLow "<<dLow<<endl;
	cout<<"dHigh "<<dHigh<<endl;

	delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsOut];
	SetBICParametersNormal(sParam,  vMean,  mWTildeOut, iNumOfObsOut,   iNumOfVar, iNumOfKnots);

	if(dHigh<dLow)
	{
		cout<<"!!!!!!!! Error dHigh < dLow !!!!!!!!!"<<endl;
	}
	double * vPredMeanOut=new double[iNumOfObsOut];
	double*  mPredDistOut=new double[iNumOfObsOut*(int(dHigh-dLow)+1)];
	BootstrapFilterOrderedNormalForecast(sParam, vDataOut,  dLow, dHigh,   iNumOfObsOut,  iNumOfParticles,vPredMeanOut,  mPredDistOut);

	string sPredDistFileOut="OutPredDist.csv";
	WriteOutDoubleArrayOut(mPredDistOut, iNumOfObsOut, ((dHigh-dLow)+1), sPrefix+sPredDistFileOut);
	string sPredMeanFileOut="OutPredMean.csv";
	WriteOutDoubleArrayOut(vPredMeanOut, iNumOfObsOut, 1, sPrefix+sPredMeanFileOut);


	cout << "Forecast Out done!"<<endl;



	delete [] vPredMeanIn;
	delete [] mPredDistIn;
	delete [] vPredMeanOut;
	delete [] mPredDistOut;
	delete [] vKnots;

	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;

	delete [] sParam.vS;
}

void OutOfSample::ForecastOrdT(int iNumOfParticles)
{

	/* Calculate W matrix for the seasonality*/
	int iNumOfKnots=3;
	double * vKnots= new double[iNumOfKnots];


//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=5400; /* 11 */
//	vKnots[2]=10800; /* 12.30 */
//	vKnots[3]=18000;/* 14.30 */
//	vKnots[iNumOfKnots-1]=23100; /* 15.55 */

	vKnots[0]=0;
	vKnots[1]=10800;
	vKnots[iNumOfKnots-1]=23400;

	gsl_matrix * mWTildeIn = gsl_matrix_alloc (iNumOfObsIn, iNumOfKnots-1);
	gsl_matrix * mWTildeOut = gsl_matrix_alloc (iNumOfObsOut, iNumOfKnots-1);
	CalculateSplineMatricesOut(iNumOfObsIn,  vTimesIn, iNumOfKnots,  vKnots,  mWTildeIn );
	CalculateSplineMatricesOut(iNumOfObsOut,  vTimesOut, iNumOfKnots,  vKnots,  mWTildeOut );

	/* Parameters */
	AllParam sParam;
	sParam.dMu=new double;
	sParam.dPhi=new double;
	sParam.dSigma2=new double;
	sParam.dGamma=new double;
	sParam.dNu=new double;
	sParam.vS=new double [iNumOfObsIn];

	cout<<"Param set start "<<endl;
	SetBICParametersT(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
	cout<<"Param set"<<endl;
	double dLow=0;
	double dHigh=0;

	for(int i=0; i<iNumOfObsOut;i++)
	{
		if(vDataIn[i]<dLow)
		{
			dLow=vDataIn[i];
		}
		if(vDataIn[i]>dHigh)
		{
			dHigh=vDataIn[i];
		}
	}

	cout<<"dLow "<<dLow<<endl;
	cout<<"dHigh "<<dHigh<<endl;

	if(dHigh<dLow)
	{
		cout<<"!!!!!!!! Error dHigh < dLow !!!!!!!!!"<<endl;
	}
	double * vPredMeanIn=new double[iNumOfObsIn];
	double*  mPredDistIn=new double[iNumOfObsIn*(int(dHigh-dLow)+1)];
	BootstrapFilterOrderedTForecast(sParam, vDataIn,  dLow, dHigh,   iNumOfObsIn,  iNumOfParticles,vPredMeanIn,  mPredDistIn);

	string sPredDistFile="InPredDist.csv";
	WriteOutDoubleArrayOut(mPredDistIn, iNumOfObsIn, ((dHigh-dLow)+1), sPrefix+sPredDistFile);
	string sPredMeanFile="InPredMean.csv";
	WriteOutDoubleArrayOut(vPredMeanIn, iNumOfObsIn, 1, sPrefix+sPredMeanFile);


	cout << "Forecast In done!"<<endl;
	dLow=0;
	dHigh=0;
	for(int i=0; i<iNumOfObsOut;i++)
	{
		if(vDataOut[i]<dLow)
		{
			dLow=vDataOut[i];
		}
		if(vDataOut[i]>dHigh)
		{
			dHigh=vDataOut[i];
		}
	}

	cout<<"dLow "<<dLow<<endl;
	cout<<"dHigh "<<dHigh<<endl;
	delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsOut];
	SetBICParametersT(sParam,  vMean,  mWTildeOut, iNumOfObsOut,   iNumOfVar, iNumOfKnots);

	if(dHigh<dLow)
	{
		cout<<"!!!!!!!! Error dHigh < dLow !!!!!!!!!"<<endl;
	}
	double * vPredMeanOut=new double[iNumOfObsOut];
	double*  mPredDistOut=new double[iNumOfObsOut*(int(dHigh-dLow)+1)];
	BootstrapFilterOrderedTForecast(sParam, vDataOut,  dLow, dHigh,   iNumOfObsOut,  iNumOfParticles,vPredMeanOut,  mPredDistOut);

	string sPredDistFileOut="OutPredDist.csv";
	WriteOutDoubleArrayOut(mPredDistOut, iNumOfObsOut, ((dHigh-dLow)+1), sPrefix+sPredDistFileOut);
	string sPredMeanFileOut="OutPredMean.csv";
	WriteOutDoubleArrayOut(vPredMeanOut, iNumOfObsOut, 1, sPrefix+sPredMeanFileOut);


	cout << "Forecast Out done!"<<endl;



	delete [] vPredMeanIn;
	delete [] mPredDistIn;
	delete [] vPredMeanOut;
	delete [] mPredDistOut;
	delete [] vKnots;

	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;
	delete sParam.dNu;
	delete [] sParam.vS;
}

OutOfSample::~OutOfSample()
{

	delete [] vDataIn;
	delete [] vTimesIn;
	delete [] vDataOut;
	delete [] vTimesOut;
	delete [] vPriceIn;
	delete [] vLogReturnIn;
	delete [] vPriceOut;
	delete [] vLogReturnOut;
	delete [] mParam;
	delete [] vMean;
	cout << "MCMC destructor has been called " << endl;
}

//void OutOfSample::ForecastLogNormal(int iNumOfParticles)
//{
//
//	/* Calculate W matrix for the seasonality*/
//	int iNumOfKnots=3;
//	double * vKnots= new double[iNumOfKnots];
//
//
////	vKnots[0]=300; /* 9.35 */
////	vKnots[1]=5400; /* 11 */
////	vKnots[2]=10800; /* 12.30 */
////	vKnots[3]=18000;/* 14.30 */
////	vKnots[iNumOfKnots-1]=23100; /* 15.55 */
//
//	vKnots[0]=0;
//	vKnots[1]=10800;
//	vKnots[iNumOfKnots-1]=23400;
//
//	gsl_matrix * mWTildeIn = gsl_matrix_alloc (iNumOfObsIn, iNumOfKnots-1);
//	gsl_matrix * mWTildeOut = gsl_matrix_alloc (iNumOfObsOut, iNumOfKnots-1);
//	CalculateSplineMatricesOut(iNumOfObsIn,  vTimesIn, iNumOfKnots,  vKnots,  mWTildeIn );
//	CalculateSplineMatricesOut(iNumOfObsOut,  vTimesOut, iNumOfKnots,  vKnots,  mWTildeOut );
//
//	/* Parameters */
//	AllParam sParam;
//	sParam.dMu=new double;
//	sParam.dPhi=new double;
//	sParam.dSigma2=new double;
//	sParam.dGamma=new double;
//	sParam.vS=new double [iNumOfObsIn];
//
//	SetBICParametersSkellam(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
//	double dRange=1;
//	double * vPredMean=new double[iNumOfObsIn];
//	double*  mPredDist=new double[iNumOfObsIn*(2*int(dRange)+1)];
//	BootstrapFilterLogNormal(sParam, vLogReturnIn,  dRange, vPriceIn,  iNumOfObsIn,  iNumOfParticles,vPredMean,  mPredDist);
//
//	string sPredDistFile="PredDist.csv";
//	WriteOutDoubleArrayOut(mPredDist, iNumOfObsIn, (2*int(dRange)+1), sPrefix+sPredDistFile);
//	string sPredMeanFile="PredMean.csv";
//	WriteOutDoubleArrayOut(vPredMean, iNumOfObsIn, 1, sPrefix+sPredMeanFile);
//
//
//	cout << "BIC In done!"<<endl;
//
//
//
//
//
//	delete [] vPredMean;
//	delete [] mPredDist;
//	delete [] vKnots;
//
//	delete sParam.dMu;
//	delete sParam.dPhi;
//	delete sParam.dSigma2;
//	delete sParam.dGamma;
//
//	delete [] sParam.vS;
//}
