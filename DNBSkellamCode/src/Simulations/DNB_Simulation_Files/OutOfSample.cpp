/*
 * OutOfSample.cpp
 *
 * The methods of the OutOfSample class
 * 
 * First created by István Barra on Aug 19, 2014
 * Modified by Agnieszka Borowska on Jan 23, 2018
 * 
 */

#include "OutOfSample.h"
#include "basic_functions2.h"
#include <time.h>

OutOfSample::OutOfSample(string sFileParam,unsigned int iBurnIn, unsigned int iThin, string sType, string sFileDataIn,string sFileDataOut )
{
	/* Import parameter draws and thin */

	CreatPrefixOut(sFileDataIn,sType, &sPrefix );
	cout<< "prefix done" <<endl;
	ImportMatrix(sPrefix+sFileParam,  iBurnIn, iThin, &mParam, &vMean, &iNumOfParam,&iNumOfVar);
	cout<< "Check" <<endl;
	cout<< sPrefix+sFileParam<<endl;
	/* Import data */
	ImportData(sFileDataIn,&vDataIn, &vTimesIn,&vPriceIn, &vLogReturnIn, &iNumOfObsIn);
	/* Import data */
	ImportData(sFileDataOut,&vDataOut, &vTimesOut,&vPriceOut, &vLogReturnOut, &iNumOfObsOut);
//	for(int i=0; i<iNumOfObsIn;i++)
//	{
//		cout<<"Data In " <<vDataIn[i]<<endl;
//	}
//	for(int i=0; i<iNumOfObsIn;i++)
//	{
//		cout<<"Data Out " <<vDataIn[i]<<endl;
//	}
}


void OutOfSample::ForecastSkellam(int iNumOfParticles)
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

	SetBICParametersSkellam(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
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
	double * vPredMean=new double[iNumOfObsIn];
	double*  mPredDist=new double[iNumOfObsIn*(int(dHigh-dLow)+1)];
	OneStepAheadSkellam(sParam, vDataIn,  dLow, dHigh, vPriceIn,  iNumOfObsIn,  iNumOfParticles,vPredMean,  mPredDist);

	string sPredDistFile="InPredDist.csv";
	WriteOutDoubleArrayOut(mPredDist, iNumOfObsIn, ((dHigh-dLow)+1), sPrefix+sPredDistFile);
	string sPredMeanFile="InPredMean.csv";
	WriteOutDoubleArrayOut(vPredMean, iNumOfObsIn, 1, sPrefix+sPredMeanFile);


	cout << "Forecast In done!"<<endl;

	delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsOut];

	SetBICParametersSkellam(sParam,  vMean,  mWTildeOut , iNumOfObsOut,   iNumOfVar, iNumOfKnots);
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
	if(dHigh<dLow)
	{
		cout<<"!!!!!!!! Error dHigh < dLow !!!!!!!!!"<<endl;
	}

	delete [] vPredMean;
	delete [] mPredDist;
	vPredMean=new double[iNumOfObsOut];
	mPredDist=new double[iNumOfObsOut*(int(dHigh-dLow)+1)];
	OneStepAheadSkellam(sParam, vDataOut,  dLow, dHigh, vPriceOut,  iNumOfObsOut,  iNumOfParticles,vPredMean,  mPredDist);

	string sPredDistFile2="OutPredDist.csv";
	WriteOutDoubleArrayOut(mPredDist, iNumOfObsOut, ((dHigh-dLow)+1), sPrefix+sPredDistFile2);
	string sPredMeanFile2="OutPredMean.csv";
	WriteOutDoubleArrayOut(vPredMean, iNumOfObsOut, 1, sPrefix+sPredMeanFile2);


	cout << "Forecast Out done!"<<endl;


	delete [] vPredMean;
	delete [] mPredDist;
	delete [] vKnots;

	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;

	delete [] sParam.vS;
}


void OutOfSample::ForecastDNB(int iNumOfParticles)
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

	SetBICParametersDNB(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
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

	double * vPredMean=new double[iNumOfObsIn];
	double*  mPredDist=new double[iNumOfObsIn*(int(dHigh-dLow)+1)];
	OneStepAheadDNB(sParam, vDataIn,  dLow, dHigh, vPriceIn,  iNumOfObsIn,  iNumOfParticles,vPredMean,  mPredDist);

	string sPredDistFile="InPredDist.csv";
	WriteOutDoubleArrayOut(mPredDist, iNumOfObsIn, ((dHigh-dLow)+1), sPrefix+sPredDistFile);
	string sPredMeanFile="InPredMean.csv";
	WriteOutDoubleArrayOut(vPredMean, iNumOfObsIn, 1, sPrefix+sPredMeanFile);


	cout << "Forecast In done!"<<endl;

	delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsOut];

	SetBICParametersDNB(sParam,  vMean,  mWTildeOut , iNumOfObsOut,   iNumOfVar, iNumOfKnots);
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
	if(dHigh<dLow)
	{
		cout<<"!!!!!!!! Error dHigh < dLow !!!!!!!!!"<<endl;
	}

	delete [] vPredMean;
	delete [] mPredDist;
	vPredMean=new double[iNumOfObsOut];
	mPredDist=new double[iNumOfObsOut*(int(dHigh-dLow)+1)];
	OneStepAheadDNB(sParam, vDataOut,  dLow, dHigh, vPriceOut,  iNumOfObsOut,  iNumOfParticles,vPredMean,  mPredDist);

	string sPredDistFile2="OutPredDist.csv";
	WriteOutDoubleArrayOut(mPredDist, iNumOfObsOut, ((dHigh-dLow)+1), sPrefix+sPredDistFile2);
	string sPredMeanFile2="OutPredMean.csv";
	WriteOutDoubleArrayOut(vPredMean, iNumOfObsOut, 1, sPrefix+sPredMeanFile2);


	cout << "Forecast Out done!"<<endl;


	delete [] vPredMean;
	delete [] mPredDist;
	delete [] vKnots;

	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;

	delete [] sParam.vS;
}

void OutOfSample::GoodnessOfFitSkellam(int iNumOfParticles)
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
		SetParametersSkellam( sParam, mParam, mWTildeOut ,iNumOfObsOut, iNumOfVar,  iNumOfKnots, i);

		/* Run particle filter */
		time_t start,end;
		time (&start);
		//BootstrapFilterSkellamAdapt( sParam,  vData, iNumOfObs, iNumOfParticles,  vTempL);
		BootstrapFilterSkellam( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
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


	string sLLFile="vLogLikeSkellam.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sLLFile);
	cout << "Predictive LL done!"<<endl;

	/* BIC Out of sample */
	SetBICParametersSkellam(sParam,  vMean,  mWTildeOut , iNumOfObsOut,   iNumOfVar, iNumOfKnots);
	BootstrapFilterSkellam( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
	for(int j=0;j<iNumOfObsOut; j++)
	{
			vLogL[j]=log(vTempL[j]);
	}

	string sOutBICFile="vOutSkellamBIC.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sOutBICFile);
	cout << "BIC Out done!"<<endl;

	/* BIC In sample */

	double * vTempLIn=new double[iNumOfObsIn];
	double * vLogLIn=new double[iNumOfObsIn];

	delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsIn];

	SetBICParametersSkellam(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
	BootstrapFilterSkellam( sParam,  vDataIn, iNumOfObsIn, iNumOfParticles,  vTempLIn);
	for(int j=0;j<iNumOfObsIn; j++)
	{
		vLogLIn[j]=log(vTempLIn[j]);
	}

	string sInBICFile="vInSkellamBIC.csv";
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

void OutOfSample::GoodnessOfFitDNB( int iNumOfParticles)
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
		SetParametersDNB( sParam, mParam, mWTildeOut ,iNumOfObsOut, iNumOfVar,  iNumOfKnots, i);

		/* Run particle filter */
		time_t start,end;
		time (&start);
		//BootstrapFilterSkellamAdapt( sParam,  vData, iNumOfObs, iNumOfParticles,  vTempL);
		BootstrapFilterDNB( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
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


	string sLLFile="vLogLikeDNB.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sLLFile);
	cout << "Predictive LL done!"<<endl;

	/* BIC Out of sample */
	SetBICParametersDNB(sParam,  vMean,  mWTildeOut , iNumOfObsOut,   iNumOfVar, iNumOfKnots);
	BootstrapFilterDNB( sParam,  vDataOut, iNumOfObsOut, iNumOfParticles,  vTempL);
	for(int j=0;j<iNumOfObsOut; j++)
	{
			vLogL[j]=log(vTempL[j]);
	}

	string sOutBICFile="vOutDNBBIC.csv";
	WriteOutDoubleArrayOut(vLogL, iNumOfObsOut, 1, sPrefix+sOutBICFile);
	cout << "BIC Out done!"<<endl;

	/* BIC In sample */

	double * vTempLIn=new double[iNumOfObsIn];
	double * vLogLIn=new double[iNumOfObsIn];

	//delete [] sParam.vS;
	sParam.vS=new double [iNumOfObsIn];

	SetBICParametersDNB(sParam,  vMean,  mWTildeIn , iNumOfObsIn,   iNumOfVar, iNumOfKnots);
	BootstrapFilterDNB( sParam,  vDataIn, iNumOfObsIn, iNumOfParticles,  vTempLIn);
	for(int j=0;j<iNumOfObsIn; j++)
	{
		vLogLIn[j]=log(vTempLIn[j]);
	}

	string sInBICFile="vInDNBBIC.csv";
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

