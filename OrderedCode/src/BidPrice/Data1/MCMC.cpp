/*
 * MCMC.cpp
 *
 *  Created on: Mar 17, 2014
 *      Author: istvan
 */

#include "MCMC.h"
#include "struct.h"
//#include "cuda_functions.h"
#include "basic_functions.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits>


MCMC::MCMC()
{

}

MCMC::MCMC(vector<double> vParam)
{

	if(vParam.size()==4 ||  vParam.size()==5 )
	{

		dMuTrue=vParam[0];
		dPhiTrue=vParam[1];
		dSigma2True=vParam[2];
		dGammaTrue=vParam[3];


		if(vParam.size()==5){
			dNuTrue=vParam[4];
		}
	}
	else
	{
		cout <<  "Number of parameters should be either 4 or 5 " << endl ;
	}

}

MCMC::~MCMC()
{
	delete [] vData;
	delete [] vTimes;
//	delete [] vSeasonTrue;
//	delete [] vXTrue;
//	delete [] vYStarTrue;
//	delete [] vLambdaTrue;
//	delete [] vNTrue;
//	delete [] vTau1True;
//	delete [] vTau2True;
	cout << "MCMC destructor has been called " << endl;
}

/*
 * Simulate  from dynamic zero inflated Skellam model
 *
 * iNumOfSim	number of simulated data points
 *
 * Note that the algorithm uses :
 * dMuTrue;
 * dPhiTrue;
 * Sigma2True;
 * dGammaTrue;
 *
 */
void MCMC::SimulateOrderedNormalData( int iNumOfSim)
{
	const gsl_rng_type * T;
	gsl_rng * r;

	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, 1);

	vData = new double[iNumOfSim];
	vXTrue = new double[iNumOfSim];
	vYStarTrue = new double[iNumOfSim];

	vSeasonTrue = new double[iNumOfSim];
	vTimes = new double[iNumOfSim];
	iNumOfObs=iNumOfSim;

	cout<<"dMuTrue "<<dMuTrue<<endl;
	cout<<"dPhiTrue "<<dPhiTrue<<endl;
	cout<<"dSigma2True "<<dSigma2True<<endl;
	cout<<"dGammaTrue "<<dGammaTrue<<endl;

	double dStateSigma=sqrt(dSigma2True);

	int iNumOfSimPerDay=2000;
	SimulateSeasonal(iNumOfSimPerDay, iNumOfSim/iNumOfSimPerDay,  vTimes,  vSeasonTrue);

	for(int i=0; i<iNumOfSim; i++){
		if(i==0)
		{

			/* Log volatility */
			vXTrue[i]= gsl_ran_gaussian(r,sqrt(dSigma2True/(1-dPhiTrue*dPhiTrue)) );

		}
		else
		{
			vXTrue [i]= dPhiTrue * vXTrue[i-1]+ gsl_ran_gaussian(r,dStateSigma) ;
		}
//		vXTrue [i]=0;
		/* Conditional Skellam return */
		double dVol= exp((dMuTrue +vSeasonTrue[i] + vXTrue[i])/2);
		double dNormal=gsl_ran_gaussian ( r, dVol);


		double dU=gsl_rng_uniform(r);
		if(dU<=dGammaTrue)
		{
//			vYStarTrue[i]=sqrt(-1);
			vYStarTrue[i]=dNormal;
			vData[i]=0;
		}
		else
		{
			vYStarTrue[i]=dNormal;
			vData[i]=round(dNormal);
			if(vData[i]==0)
			{
				vData[i]=0;
			}
		}

		cout << "Data at "<< i<<" is " << vData[i]<<  endl;
	}

	gsl_rng_free(r);
	return;
}

void MCMC::SimulateOrderedTData( int iNumOfSim)
{
	const gsl_rng_type * T;
	gsl_rng * r;

	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, 1);

	vData = new double[iNumOfSim];
	vXTrue = new double[iNumOfSim];
	vYStarTrue = new double[iNumOfSim];
	vLambdaTrue = new double[iNumOfSim];

	vSeasonTrue = new double[iNumOfSim];
	vTimes = new double[iNumOfSim];
	iNumOfObs=iNumOfSim;

	cout<<"dMuTrue "<<dMuTrue<<endl;
	cout<<"dPhiTrue "<<dPhiTrue<<endl;
	cout<<"dSigma2True "<<dSigma2True<<endl;
	cout<<"dGammaTrue "<<dGammaTrue<<endl;
	cout<<"dNuTrue "<<dNuTrue<<endl;

	double dStateSigma=sqrt(dSigma2True);

	int iNumOfSimPerDay=2000;
	SimulateSeasonal(iNumOfSimPerDay, iNumOfSim/iNumOfSimPerDay,  vTimes,  vSeasonTrue);



    //double dRoundTo=100000000;
	double dRoundTo=1;

	for(int i=0; i<iNumOfSim; i++){
		if(i==0)
		{

			/* Log volatility */
			vXTrue[i]= gsl_ran_gaussian(r,sqrt(dSigma2True/(1-dPhiTrue*dPhiTrue)) );

		}
		else
		{
			vXTrue [i]= dPhiTrue * vXTrue[i-1]+ gsl_ran_gaussian(r,dStateSigma) ;
		}
//		vXTrue [i]=0;
		/* Conditional Skellam return */
		double dVol= exp((dMuTrue +vSeasonTrue[i] + vXTrue[i])/2);
		double dNormal=gsl_ran_gaussian ( r, dVol);
		double dScale=sqrt(1.0/gsl_ran_gamma(r,dNuTrue/2,2/dNuTrue ));


		double dU=gsl_rng_uniform(r);
		if(dU<=dGammaTrue)
		{
//			vYStarTrue[i]=sqrt(-1);
			vYStarTrue[i]=dScale*dNormal;
			vLambdaTrue[i]=dScale*dScale;
			vData[i]=0;
		}
		else
		{
			vYStarTrue[i]=dScale*dNormal;
			vLambdaTrue[i]=dScale*dScale;
			vData[i]=round(dScale*dNormal*dRoundTo)/dRoundTo;
			if(vData[i]==0)
			{
				vData[i]=0;
			}
		}

		cout << "Y star at "<< i<<" is " << vYStarTrue[i]<<  endl;
		cout << "Data at "<< i<<" is " << vData[i]<<  endl;
	}

	gsl_rng_free(r);
	return;
}

void MCMC::ImportData(string sFile)
{
	vector<vector<double> > mData;

	sInput=sFile;

	fstream ifsDataFile;
	ifsDataFile.open(sFile.c_str(), ios::in);
	if(!ifsDataFile.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	string line;
	iNumOfObs=0;
	while ( getline (ifsDataFile,line) )
	{
		vector<double>vRow;
	//	cout << line << '\n';
		stringstream ssLine(line);
		string var;
		char* end;
		while(getline(ssLine, var,',') )
		{
			vRow.push_back(strtod( var.c_str() ,&end) );
		}
		mData.push_back(vRow);

		iNumOfObs=iNumOfObs+1;
		vRow.clear();
	}

	ifsDataFile.close();

	vTimes = new double[iNumOfObs];
	vData = new double[iNumOfObs];
	
	int k=0;
	for(vector<vector<double> >::iterator i=mData.begin(); i != mData.end(); ++i)
	{
		vTimes[k]= (*i)[0];
		vData[k]= (*i)[1];
//		cout<<"Import Time at "<< k << " is "<< vTimes[k]<<endl;
//		cout<<"Import Data at "<< k << " is "<< vData[k]<<endl;
		k=k+1;
	}

	cout <<"The number of observations is: "<<iNumOfObs<<endl;
}


void MCMC::ImportData2(string sPath, string sFile)
{
	vector<vector<double> > mData;

	sInput=sFile;
	
	string sAddress;
	sAddress = sPath + sFile; 

	fstream ifsDataFile;
	// ifsDataFile.open(sFile.c_str(), ios::in);
	ifsDataFile.open(sAddress.c_str(), ios::in);
	if(!ifsDataFile.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	string line;
	iNumOfObs=0;
	while ( getline (ifsDataFile,line) )
	{
		vector<double>vRow;
	//	cout << line << '\n';
		stringstream ssLine(line);
		string var;
		char* end;
		while(getline(ssLine, var,',') )
		{
			vRow.push_back(strtod( var.c_str() ,&end) );
		}
		mData.push_back(vRow);

		iNumOfObs=iNumOfObs+1;
		vRow.clear();
	}
	ifsDataFile.close();

	vTimes = new double[iNumOfObs];
	vData = new double[iNumOfObs];
	
	int k=0;
	for(vector<vector<double> >::iterator i=mData.begin(); i != mData.end(); ++i)
	{
		vTimes[k]= (*i)[0];
		vData[k]= (*i)[1];
		cout<<"vData at "<<k<<" is "<<vData[k]<<endl;
		k=k+1;
	}
	cout <<"The number of observations is: "<<iNumOfObs<<endl;
}


void MCMC::ExportSimulatedData()
{
	fstream ofsData;
	ofsData.open("SimulatedData.csv", ios::out);
	if(!ofsData.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	for(int i=0; i< iNumOfObs; ++i){
		ofsData << vTimes[i]<<" , "<< vData[i]<<" , "<<1<<" , "<<1<< endl;
	}
	ofsData.close();

	fstream ofsLogInt;
	ofsLogInt.open("SimulatedX.csv", ios::out);
	if(!ofsLogInt.is_open()){
		cerr << "Cannot open the file" << endl;
		return;
	}
	for(int i=0; i< iNumOfObs; ++i){
		ofsLogInt << vXTrue[i] << endl;
	}
	ofsLogInt.close();

	return;
}

void MCMC::ExportEstimationResults()
{
	string sPrefix;
	CreatPrefix(sInput,sType, &sPrefix );


	fstream ofsData;
	ofsData.open((sPrefix+"EstimationResults.csv").c_str(), ios::out);
	if(!ofsData.is_open()){
		cerr << " ExportEstimationResults Error: Cannot open the file" << endl;
		return;
	}
	for(vector<vector<double> >::iterator i=mChain.begin(); i != mChain.end(); ++i){

		for(vector<double>::iterator j= i->begin(); j != i -> end(); ++j){
			if(j==i->end()-1)
			{
				ofsData << *j  << endl;
			}
			else
			{
				ofsData << *j  << ",";
			}
		}
	}
	ofsData.close();
}

void MCMC::EstimateNormal( int iNumOfIter,int iBurnIn)
{
	sType="OrdNormal";
	CreatPrefix(sInput,sType, &sPrefix );
	cout<< "prefix done" <<endl;
	/* Initialize parameters and latent variables */
	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	int iSeedNo = 2;
	gsl_rng_set(gsl_random_num, iSeedNo);



	double* dCovar = new double;
	double* dSum=new double;
	double* mCovar = new double[4];
	double* mSum=new double[2];
	mCovar[0]=0;
	mCovar[1]=0;
	mCovar[2]=0;
	mCovar[3]=0;
	mSum[0]=0;
	mSum[1]=0;
	dCovar[0]=0;
	dSum[0]=0;

	/* Initialize Parameters*/
	AllParam sParam;
	sParam.dMu=new double;
	sParam.iNumOfObs=new int;
	sParam.dPhi=new double;
	sParam.dSigma2=new double;
	sParam.dGamma=new double;
	sParam.dPriorGammaA=new double;
	sParam.dPriorGammaB=new double;
	sParam.dPriorMuMean=new double;
	sParam.dPriorMuSigma2=new double;
	sParam.dPriorPhiA=new double;
	sParam.dPriorPhiB=new double;
	sParam.dPriorSigmaA=new double;
	sParam.dPriorSigmaB=new double;
	sParam.dPriorSeasonalMean=new double;
	sParam.dPriorSeasonalVar=new double;

	sParam.vX=new double [iNumOfObs];
	sParam.vS=new double [iNumOfObs];
	sParam.mAuxY=new double[iNumOfObs];
	sParam.mAuxH=new double[iNumOfObs];
	sParam.vYStar=new double[iNumOfObs];
	sParam.vLambda=new double[iNumOfObs];


	sParam.iNumOfObs[0]=iNumOfObs;

	sParam.dMu[0]= 3.4; //0;
	sParam.dPhi[0]= 0.97; //0.95;
	sParam.dSigma2[0]= 0.065; //0.01;
	sParam.dGamma[0]=0.45; //0.3; //0;




	sParam.dPriorGammaA[0]=1.7;
	sParam.dPriorGammaB[0]=10;
	sParam.dPriorMuMean[0]=0;
	sParam.dPriorMuSigma2[0]=10; //1;
	sParam.dPriorPhiA[0]=20;
	sParam.dPriorPhiB[0]=1.5;
	sParam.dPriorSigmaA[0]=2.5;
	sParam.dPriorSigmaB[0]=1/0.025;
	sParam.dPriorSeasonalMean[0]=0;
	sParam.dPriorSeasonalVar[0]=1;


//	for(int i=0; i<iNumOfObs; i++)
//	{
//
//		sParam.vYStar[i]=vYStarTrue[i];
//		sParam.vS[i]=vSeasonTrue[i];
//		sParam.vX[i]=vXTrue[i];
//
//	}


	/*Initializing the iteration*/


	time_t start,end;
	time (&start);



	double * vLogIntEst=new double[iNumOfObs];
	double * vXEst=new double[iNumOfObs];
	double * vSeasonEst=new double[iNumOfObs];
	double * vYStarEst=new double[iNumOfObs];
	double * vLambdaEst=new double[iNumOfObs];
	double * vVolEst=new double[iNumOfObs];
	double * vX2Est=new double[iNumOfObs];
	double * vSeason2Est=new double[iNumOfObs];
	double * vVol2Est=new double[iNumOfObs];
	double * vXVolEst=new double[iNumOfObs];
	double * vSeasonVolEst=new double[iNumOfObs];
	double * vVolVolEst=new double[iNumOfObs];

	for(int j=0; j<iNumOfObs;j++)
	{
		vXEst[j]=0;
		vSeasonEst[j]=0;
		vVolEst[j]=0;
		vX2Est[j]=0;
		vSeason2Est[j]=0;
		vVol2Est[j]=0;
		vLogIntEst[j]=0;
		vYStarEst[j]=0;
		vLambdaEst[j]=0;

	}

//
//	double * mXDraws=new double[iNumOfObs*iNumOfIter];
//	double * mNDraws=new double[iNumOfObs*iNumOfIter];
//	double * mSeasonDraws=new double[iNumOfObs*iNumOfIter];
//	double * mLogIntDraws=new double[iNumOfObs*iNumOfIter];
//	double * mTau1Draws=new double[iNumOfObs*iNumOfIter];
//	double * mTau2Draws=new double[iNumOfObs*iNumOfIter];

//	string svYTrueFile="mYTrue.csv";
//	WriteOutDoubleArray( vData , iNumOfObs, 1, svYTrueFile);


	int iNumOfKnots=3;
	double * vKnots= new double[iNumOfKnots];
//	vKnots[0]=0; /* 9.30 */
//	vKnots[1]=5400; /* 11 */
//	vKnots[2]=10800; /* 12.30 */
//	vKnots[3]=18000;/* 14.30 */
//	vKnots[iNumOfKnots-1]=23400; /* 16.00 */

//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=5400; /* 11 */
//	vKnots[2]=10800; /* 12.30 */
//	vKnots[3]=18000;/* 14.30 */
//	vKnots[iNumOfKnots-1]=23100; /* 15.55 */

//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=1800;
//	vKnots[2]=3800;
//	vKnots[3]=5400; /* 11 */
//	vKnots[4]=6700;
//	vKnots[5]=8500;
//	vKnots[6]=10800; /* 12.30 */
//	vKnots[7]=14000;
//	vKnots[8]=18000;/* 14.30 */
//	vKnots[9]=20000;
//	vKnots[iNumOfKnots-1]=23100; /* 15.55 */

//	vKnots[0]=300; /* 9.35 */
//	vKnots[1]=23100; /* 15.55 */

	vKnots[0]=0; /* 9.30 */
	vKnots[1]=10800;  /* 12.30 */
	vKnots[2]=23400; /* 16.00 */
	cout<<"iNumOfKnots "<<iNumOfKnots<<endl;
//	vKnots[0]=0; /* 9.30 */
//	vKnots[1]=10800;
//	vKnots[2]=23400;  /* 16.00 */

	sParam.vBeta=new double [iNumOfKnots-1];
	for(int k=0; k<iNumOfKnots-1; k++)
	{
		sParam.vBeta[k]=0;
	}

	gsl_matrix * mWTilde = gsl_matrix_alloc (iNumOfObs, iNumOfKnots-1);
	CalculateSplineMatrices(iNumOfObs,  vTimes, iNumOfKnots,  vKnots,  mWTilde );

	double * vInitialLogVol=new double[iNumOfObs];
	CalculateInitialLogVol(vData,  vTimes,  iNumOfObs, 50, vInitialLogVol);
	CalculateInitialState(sParam,  vInitialLogVol,  iNumOfObs, iNumOfKnots,   mWTilde  );

	sParam.iNumOfKnots=new  int;
	sParam.iNumOfKnots[0]=iNumOfKnots;

	sParam.W=new gsl_matrix;
	sParam.W=mWTilde;

	string svInitLogIntFile="vInitialLogVol.csv";
	WriteOutDoubleArray( vInitialLogVol, iNumOfObs, 1, svInitLogIntFile);
	string svInitSFile="vInitS.csv";
	WriteOutDoubleArray( sParam.vS, iNumOfObs, 1, svInitSFile);
	string svInitXFile="vInitX.csv";
	WriteOutDoubleArray( sParam.vX, iNumOfObs, 1, svInitXFile);


	double * mMarkerXLow=new double [iNumOfObs*5];
	double * mMarkerSLow=new double [iNumOfObs*5];
	double * mMarkerVolLow=new double [iNumOfObs*5];
	double * mMarkerXHigh=new double [iNumOfObs*5];
	double * mMarkerSHigh=new double [iNumOfObs*5];
	double * mMarkerVolHigh=new double [iNumOfObs*5];

	double * mNXLow=new double [iNumOfObs*5];
	double * mNSLow=new double [iNumOfObs*5];
	double * mNVolLow=new double [iNumOfObs*5];
	double * mNXHigh=new double [iNumOfObs*5];
	double * mNSHigh=new double [iNumOfObs*5];
	double * mNVolHigh=new double [iNumOfObs*5];

	double * mNDXLow=new double [iNumOfObs*5];
	double * mNDSLow=new double [iNumOfObs*5];
	double * mNDVolLow=new double [iNumOfObs*5];
	double * mNDXHigh=new double [iNumOfObs*5];
	double * mNDSHigh=new double [iNumOfObs*5];
	double * mNDVolHigh=new double [iNumOfObs*5];

	double * vVol=new double [iNumOfObs];

	cout<<"start"<<endl;
	for( int i=0; i<iNumOfIter; i++)
	{
		cout <<"==== Iteration: " << i << " ===="<< endl;
		/*Draw vY Star */
		DrawYStarNormal( sParam, vData,iNumOfObs ,  gsl_random_num);
		cout << "YStar is Done!" <<endl;
		/*Draw AuxY AuxH*/
		DrawAuxYandAuxHNormal(sParam, iNumOfObs , vData, gsl_random_num);
		cout<<"Aux Done "<<endl;

		/* Draw Gamma */
		DrawGammaAdaptiveRWNormal(vData,  iNumOfObs,i, gsl_random_num, sParam,dCovar,  dSum  );
//		cout << "dGamma " << sParam.dGamma[0] << endl;



		/*Draw Phi, Sigma2*/
		DrawPhiSigmaAdaptiveRWWithSeasonal( iNumOfKnots, mWTilde,sParam, iNumOfObs,  i , gsl_random_num,  mCovar,   mSum );
		cout << "dPhi " << sParam.dPhi[0] << endl;
		cout << "dSigma " << sParam.dSigma2[0] << endl;

		/*Draw X, Mu*/

		DrawXandMuWithSeasonal(iNumOfKnots, mWTilde, sParam,iNumOfObs, gsl_random_num );
		cout << "dMu " << sParam.dMu[0] << endl;


		if(i>=iBurnIn)
		{
			for(int k=0; k<iNumOfObs; k++)
			{
				vVol[k]=sqrt((1-sParam.dGamma[0])*VarOrdNormal(exp((sParam.dMu[0]+sParam.vS[k]+sParam.vX[k])/2)));
			}

			if(i<iBurnIn+5)
			{
				for(int k=0; k<iNumOfObs; k++)
				{
					mMarkerXLow[k*5+i-iBurnIn]=sParam.vX[k];
					mMarkerSLow[k*5+i-iBurnIn]=sParam.vS[k];;
					mMarkerVolLow[k*5+i-iBurnIn]=vVol[k];

					mMarkerXHigh[k*5+i-iBurnIn]=sParam.vX[k];
					mMarkerSHigh[k*5+i-iBurnIn]=sParam.vS[k];;
					mMarkerVolHigh[k*5+i-iBurnIn]=vVol[k];
				}

			}
			else
			{
				if(i==iBurnIn+5)
				{
					InitializeQuantilePPAlgo(iNumOfObs, 0.05, mMarkerXLow,  mNXLow , mNDXLow);
					InitializeQuantilePPAlgo(iNumOfObs, 0.05, mMarkerSLow,  mNSLow , mNDSLow);
					InitializeQuantilePPAlgo(iNumOfObs, 0.05, mMarkerVolLow,  mNVolLow , mNDVolLow);
					InitializeQuantilePPAlgo(iNumOfObs, 0.95, mMarkerXHigh,  mNXHigh , mNDXHigh);
					InitializeQuantilePPAlgo(iNumOfObs, 0.95, mMarkerSHigh,  mNSHigh , mNDSHigh);
					InitializeQuantilePPAlgo(iNumOfObs, 0.95, mMarkerVolHigh,  mNVolHigh , mNDVolHigh);
				}

				 QuantilePPAlgo(iNumOfObs, 0.05,sParam.vX,  mMarkerXLow,  mNXLow , mNDXLow);
				 QuantilePPAlgo(iNumOfObs, 0.05, sParam.vS, mMarkerSLow,  mNSLow , mNDSLow);
				 QuantilePPAlgo(iNumOfObs, 0.05,vVol, mMarkerVolLow,  mNVolLow , mNDVolLow);
				 QuantilePPAlgo(iNumOfObs, 0.95,sParam.vX, mMarkerXHigh,  mNXHigh , mNDXHigh);
				 QuantilePPAlgo(iNumOfObs, 0.95,sParam.vS, mMarkerSHigh,  mNSHigh , mNDSHigh);
				 QuantilePPAlgo(iNumOfObs, 0.95,vVol, mMarkerVolHigh,  mNVolHigh , mNDVolHigh);
			}


//			SequentialStore(i-iBurnIn, iNumOfIter-iBurnIn,iNumOfObs, sParam.vX , Xfs,  sPrefix+sXDraw);
//			SequentialStore(i-iBurnIn, iNumOfIter-iBurnIn,iNumOfObs, sParam.vS , Sfs,  sPrefix+sSDraw);

			for(int j=0; j<iNumOfObs;j++)
			{
				vXEst[j]=vXEst[j]+(sParam.vX[j])/(iNumOfIter-iBurnIn);
				vSeasonEst[j]=vSeasonEst[j]+(sParam.vS[j])/(iNumOfIter-iBurnIn);
				vVolEst[j]=vVolEst[j]+vVol[j]/(iNumOfIter-iBurnIn);
				vLogIntEst[j]=vLogIntEst[j]+(sParam.dMu[0]+sParam.vS[j]+sParam.vX[j])/(iNumOfIter-iBurnIn);

				vYStarEst[j]=vYStarEst[j]+(sParam.vYStar[j])/(iNumOfIter-iBurnIn);
				vLambdaEst[j]=vLambdaEst[j]+(sParam.vLambda[j])/(iNumOfIter-iBurnIn);

				vX2Est[j]=vX2Est[j]+sParam.vX[j]*sParam.vX[j]/(iNumOfIter-iBurnIn);
				vSeason2Est[j]=vSeason2Est[j]+sParam.vS[j]*sParam.vS[j]/(iNumOfIter-iBurnIn);
				vVol2Est[j]=vVol2Est[j]+vVol[j]*vVol[j]/(iNumOfIter-iBurnIn);
			}
		}


		/*Save new draws */
		vector<double> vTempRow;

		vTempRow.push_back(sParam.dMu[0]);
		vTempRow.push_back(sParam.dPhi[0]);
		vTempRow.push_back(sParam.dSigma2[0]);
		vTempRow.push_back(sParam.dGamma[0]);
		for(int k=0; k<iNumOfKnots-1; k++)
		{
			vTempRow.push_back(sParam.vBeta[k]);
		}

		mChain.push_back(vTempRow);
	}
	time (&end);
	double dif = difftime (end,start);
	printf ("Elasped time is %.2lf seconds. \n", dif );


		FILE * pFile;
		pFile = fopen ("OrdN_time.txt","w");
		fprintf (pFile, "OrdN time = %16.4f s\n",dif);
		fclose (pFile);
		
		
	for(int j=0; j<iNumOfObs;j++)
	{
		vXVolEst[j]=vX2Est[j]-vXEst[j]*vXEst[j];
		vSeasonVolEst[j]=vSeason2Est[j]-vSeasonEst[j]*vSeasonEst[j];
		vVolVolEst[j]=vVol2Est[j]-vVolEst[j]*vVolEst[j];
	}

//	string sXTrueFile="vXTrue.csv";
//	WriteOutDoubleArray(vXTrue, iNumOfObs, 1, sXTrueFile);

	string sXEstFile="vXEst.csv";
	WriteOutDoubleArray( vXEst , iNumOfObs, 1, sPrefix+sXEstFile);
	cout<<"X Est Saved!"<<endl;

//	string sLogIntTrueFile="vLogIntTrue.csv";
//	WriteOutDoubleArray(vLogIntTrue, iNumOfObs, 1, sLogIntTrueFile);

	string sLogIntEstFile="vLogIntEst.csv";
	WriteOutDoubleArray( vLogIntEst , iNumOfObs, 1, sPrefix+sLogIntEstFile);
	cout<<"LogInt Est Saved!"<<endl;

	string sSeasonEstFile="vSeasonEst.csv";
	WriteOutDoubleArray( vSeasonEst , iNumOfObs, 1, sPrefix+sSeasonEstFile);
	cout<<"Season Est Saved!"<<endl;


	string sYStarEstFile="vYStarEst.csv";
	WriteOutDoubleArray( vYStarEst , iNumOfObs, 1, sPrefix+sYStarEstFile);
	cout<<"YStar Est Saved!"<<endl;

//	string sYStarTrueFile="vYStarTrue.csv";
//	WriteOutDoubleArray( vYStarTrue , iNumOfObs, 1, sYStarTrueFile);
//
//	string sLambdaTrueFile="vLambdaTrue.csv";
//	WriteOutDoubleArray( vLambdaTrue , iNumOfObs, 1, sLambdaTrueFile);

	string sLambdaEstFile="vLambdaEst.csv";
	WriteOutDoubleArray(vLambdaEst , iNumOfObs, 1, sPrefix+sLambdaEstFile);
	cout<<"Lambda Est Saved!"<<endl;

	string sVolEstFile="vVolEst.csv";
	WriteOutDoubleArray( vVolEst , iNumOfObs, 1, sPrefix+sVolEstFile);
	cout<<"Vol Est Saved!"<<endl;

	string sXLowFile="vXLow.csv";
	WriteOutDoubleArray(  mMarkerXLow, iNumOfObs, 5, sPrefix+sXLowFile);
	cout<<"X Low Saved!"<<endl;

	string sXHighFile="vXHigh.csv";
	WriteOutDoubleArray(  mMarkerXHigh, iNumOfObs, 5, sPrefix+sXHighFile);
	cout<<"X High Saved!"<<endl;

	string sSLowFile="vSLow.csv";
	WriteOutDoubleArray(  mMarkerSLow, iNumOfObs, 5, sPrefix+sSLowFile);
	cout<<"S Low Saved!"<<endl;

	string sSHighFile="vSHigh.csv";
	WriteOutDoubleArray(  mMarkerSHigh, iNumOfObs, 5, sPrefix+sSHighFile);
	cout<<"S High Saved!"<<endl;

	string sVolLowFile="vVolLow.csv";
	WriteOutDoubleArray(  mMarkerVolLow, iNumOfObs, 5, sPrefix+sVolLowFile);
	cout<<"Vol Low Saved!"<<endl;

	string sVolHighFile="vVolHigh.csv";
	WriteOutDoubleArray(  mMarkerVolHigh, iNumOfObs, 5, sPrefix+sVolHighFile);
	cout<<"Vol High Saved!"<<endl;

	string sXVolFile="vXVol.csv";
	WriteOutDoubleArray(  vXVolEst, iNumOfObs, 1, sPrefix+sXVolFile);
	cout<<"X Vol Saved!"<<endl;

	string sSVolFile="vSVol.csv";
	WriteOutDoubleArray(  vSeasonVolEst, iNumOfObs, 1, sPrefix+sSVolFile);
	cout<<"S Vol Saved!"<<endl;

	string sVolVolFile="vVolVol.csv";
	WriteOutDoubleArray(  vVolVolEst, iNumOfObs, 1, sPrefix+sVolVolFile);
	cout<<"Vol Vol Saved!"<<endl;

	cout <<"Write out done"<<endl;


	delete [] sParam.vX;
	delete [] sParam.vBeta;
	delete [] sParam.vS;
	delete [] sParam.mAuxY;
	delete [] sParam.mAuxH;
	delete [] sParam.vYStar;
	delete [] sParam.vLambda;

	delete [] vInitialLogVol;

	delete [] vLogIntEst;
	delete [] vXEst;
	delete [] vSeasonEst;
	delete [] vVolEst;
	delete [] vX2Est;
	delete [] vSeason2Est;
	delete [] vVol2Est;
	delete [] vYStarEst;
	delete [] vXVolEst;
	delete [] vSeasonVolEst;
	delete [] vVolVolEst;

	delete [] mMarkerXLow;
	delete [] mMarkerSLow;
	delete [] mMarkerVolLow;
	delete [] mMarkerXHigh;
	delete [] mMarkerSHigh;
	delete [] mMarkerVolHigh;

	delete [] mNXLow;
	delete [] mNSLow;
	delete [] mNVolLow;
	delete [] mNXHigh;
	delete [] mNSHigh;
	delete [] mNVolHigh;

	delete [] vVol;

	delete [] mNDXLow;
	delete [] mNDSLow;
	delete [] mNDVolLow;
	delete [] mNDXHigh;
	delete [] mNDSHigh;
	delete [] mNDVolHigh;

	delete [] vKnots;

	delete dCovar;
	delete dSum;
	delete [] mCovar;
	delete [] mSum;

	delete sParam.iNumOfObs;
	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;
	delete sParam.dPriorGammaA;
	delete sParam.dPriorGammaB;
	delete sParam.dPriorMuMean;
	delete sParam.dPriorMuSigma2;
	delete sParam.dPriorPhiA;
	delete sParam.dPriorPhiB;
	delete sParam.dPriorSigmaA;
	delete sParam.dPriorSigmaB;

	cout <<"Clear memory done"<<endl;
}

void MCMC::EstimateT( int iNumOfIter,int iBurnIn)
{
	sType="OrdT";
	CreatPrefix(sInput,sType, &sPrefix );
	cout<< "prefix done" <<endl;
	/* Initialize parameters and latent variables */
	gsl_rng_env_setup();
	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 41);



	double* dCovar = new double;
	double* dSum=new double;
	double* mCovar = new double[4];
	double* mSum=new double[2];
	mCovar[0]=0;
	mCovar[1]=0;
	mCovar[2]=0;
	mCovar[3]=0;
	mSum[0]=0;
	mSum[1]=0;
	dCovar[0]=0;
	dSum[0]=0;

	/* Initialize Parameters*/
	AllParam sParam;
	sParam.dMu=new double;
	sParam.iNumOfObs=new int;
	sParam.dPhi=new double;
	sParam.dSigma2=new double;
	sParam.dGamma=new double;
	sParam.dNu=new double;
	sParam.dPriorGammaA=new double;
	sParam.dPriorGammaB=new double;
	sParam.dPriorMuMean=new double;
	sParam.dPriorNuA=new double;
	sParam.dPriorNuB=new double;
	sParam.dPriorMuSigma2=new double;
	sParam.dPriorPhiA=new double;
	sParam.dPriorPhiB=new double;
	sParam.dPriorSigmaA=new double;
	sParam.dPriorSigmaB=new double;
	sParam.dPriorSeasonalMean=new double;
	sParam.dPriorSeasonalVar=new double;

	sParam.vX=new double [iNumOfObs];
	sParam.vS=new double [iNumOfObs];
	sParam.mAuxY=new double[iNumOfObs];
	sParam.mAuxH=new double[iNumOfObs];
	sParam.vYStar=new double[iNumOfObs];
	sParam.vLambda=new double[iNumOfObs];

	sParam.iNumOfObs[0]=iNumOfObs;

	sParam.dMu[0]=0;
	sParam.dPhi[0]=0.95;
	sParam.dSigma2[0]=0.01;
	sParam.dGamma[0]=0.3;
	sParam.dNu[0]=20;



	sParam.dPriorGammaA[0]=1.7;
	sParam.dPriorGammaB[0]=10;
	sParam.dPriorNuA[0]=9;
	sParam.dPriorNuB[0]=1.5;
	sParam.dPriorMuMean[0]=0;
	sParam.dPriorMuSigma2[0]=1;
	sParam.dPriorPhiA[0]=20;
	sParam.dPriorPhiB[0]=1.5;
	sParam.dPriorSigmaA[0]=2.5;
	sParam.dPriorSigmaB[0]=1/0.025;
	sParam.dPriorSeasonalMean[0]=0;
	sParam.dPriorSeasonalVar[0]=1;


	for(int i=0; i<iNumOfObs; i++)
	{
//
//		sParam.vYStar[i]=vYStarTrue[i];
//		sParam.vS[i]=vSeasonTrue[i];
//		sParam.vX[i]=vXTrue[i];
		sParam.vLambda[i]=1;

	}


	/*Initializing the iteration*/


	time_t start,end;
	time (&start);



	double * vLogIntEst=new double[iNumOfObs];
	double * vXEst=new double[iNumOfObs];
	double * vSeasonEst=new double[iNumOfObs];
	double * vYStarEst=new double[iNumOfObs];
	double * vLambdaEst=new double[iNumOfObs];
	double * vVolEst=new double[iNumOfObs];
	double * vX2Est=new double[iNumOfObs];
	double * vSeason2Est=new double[iNumOfObs];
	double * vVol2Est=new double[iNumOfObs];
	double * vXVolEst=new double[iNumOfObs];
	double * vSeasonVolEst=new double[iNumOfObs];
	double * vVolVolEst=new double[iNumOfObs];

	for(int j=0; j<iNumOfObs;j++)
	{


		vXEst[j]=0;
		vSeasonEst[j]=0;
		vVolEst[j]=0;
		vX2Est[j]=0;
		vSeason2Est[j]=0;
		vVol2Est[j]=0;
		vLogIntEst[j]=0;
		vYStarEst[j]=0;
		vLambdaEst[j]=0;

	}

//	int iNumOfT=10000;
//	double * vTestT=new double[iNumOfT];
//	double * vTestN=new double[iNumOfT];
//	for(int j=0;j<iNumOfT;j++)
//	{
//		vTestT[j]=DrawTruncatedT( 0.5, 10, 15.5,16.5 , gsl_random_num);
//		vTestN[j]=DrawTruncatedNormal( 0, 1, -9.5,-8.5 , gsl_random_num);
//	}
//	string svTFile="vTestT.csv";
//	WriteOutDoubleArray( vTestT, iNumOfT, 1, svTFile);
//
//	string svNFile="vTestN.csv";
//	WriteOutDoubleArray( vTestN, iNumOfT, 1, svNFile);
//
//	delete [] vTestT;
//	delete [] vTestN;
//
//	double * mXDraws=new double[iNumOfObs*iNumOfIter];
//	double * mNDraws=new double[iNumOfObs*iNumOfIter];
//	double * mSeasonDraws=new double[iNumOfObs*iNumOfIter];
//	double * mLogIntDraws=new double[iNumOfObs*iNumOfIter];
//	double * mTau1Draws=new double[iNumOfObs*iNumOfIter];
//	double * mTau2Draws=new double[iNumOfObs*iNumOfIter];

//	string svYTrueFile="mYTrue.csv";
//	WriteOutDoubleArray( vData , iNumOfObs, 1, svYTrueFile);


	int iNumOfKnots=3;
	double * vKnots= new double[iNumOfKnots];

	vKnots[0]=0; /* 9.30 */
	vKnots[1]=10800;  /* 12.30 */
	vKnots[2]=23400; /* 16.00 */
	cout<<"iNumOfKnots "<<iNumOfKnots<<endl;


	sParam.vBeta=new double [iNumOfKnots-1];
	for(int k=0; k<iNumOfKnots-1; k++)
	{
		sParam.vBeta[k]=0;
	}

	gsl_matrix * mWTilde = gsl_matrix_alloc (iNumOfObs, iNumOfKnots-1);
	CalculateSplineMatrices(iNumOfObs,  vTimes, iNumOfKnots,  vKnots,  mWTilde );

	double * vInitialLogVol=new double[iNumOfObs];
	CalculateInitialLogVol(vData,  vTimes,  iNumOfObs, 50, vInitialLogVol);
	CalculateInitialState(sParam,  vInitialLogVol,  iNumOfObs, iNumOfKnots,   mWTilde  );

	sParam.iNumOfKnots=new  int;
	sParam.iNumOfKnots[0]=iNumOfKnots;

	sParam.W=new gsl_matrix;
	sParam.W=mWTilde;

//	string svInitLogIntFile="vInitialLogVol.csv";
//	WriteOutDoubleArray( vInitialLogVol, iNumOfObs, 1, svInitLogIntFile);
//	string svInitSFile="vInitS.csv";
//	WriteOutDoubleArray( sParam.vS, iNumOfObs, 1, svInitSFile);
//	string svInitXFile="vInitX.csv";
//	WriteOutDoubleArray( sParam.vX, iNumOfObs, 1, svInitXFile);

//	int iNumOfQ=100000;
//
//	double * vTestQ=new double[iNumOfQ];
//	double * mMarker=new double[5];
//	double * mN=new double[5];
//	double * mND=new double[5];
//
//	for(int i=0; i<iNumOfQ;i++)
//	{
//		vTestQ[i]=gsl_ran_gaussian(gsl_random_num ,1);
//		if(i<5)
//		{
//			mMarker[i]=vTestQ[i];
//		}
//	}
//
//
//	InitializeQuantilePPAlgo(1, 0.05,  mMarker ,  mN ,  mND);
//	for( int i=5; i<iNumOfQ; i++)
//	{
//		QuantilePPAlgo(1, 0.05,&vTestQ[i],  mMarker,  mN , mND);
//		cout <<"Min "<<mMarker[0]<<" p "<<mMarker[1]<<" Quantile "<<mMarker[2]<<" 2p "<<mMarker[3]<<" Max "<<mMarker[4]<<endl;
//	}
//
//	cout <<"Quantile "<<mMarker[2]<<endl;
//	delete [] vTestQ;
//	delete [] mMarker;
//	delete [] mN;
//	delete [] mND;
//
//	cout<<"test end "<<endl;

	double * mMarkerXLow=new double [iNumOfObs*5];
	double * mMarkerSLow=new double [iNumOfObs*5];
	double * mMarkerVolLow=new double [iNumOfObs*5];
	double * mMarkerXHigh=new double [iNumOfObs*5];
	double * mMarkerSHigh=new double [iNumOfObs*5];
	double * mMarkerVolHigh=new double [iNumOfObs*5];

	double * mNXLow=new double [iNumOfObs*5];
	double * mNSLow=new double [iNumOfObs*5];
	double * mNVolLow=new double [iNumOfObs*5];
	double * mNXHigh=new double [iNumOfObs*5];
	double * mNSHigh=new double [iNumOfObs*5];
	double * mNVolHigh=new double [iNumOfObs*5];

	double * mNDXLow=new double [iNumOfObs*5];
	double * mNDSLow=new double [iNumOfObs*5];
	double * mNDVolLow=new double [iNumOfObs*5];
	double * mNDXHigh=new double [iNumOfObs*5];
	double * mNDSHigh=new double [iNumOfObs*5];
	double * mNDVolHigh=new double [iNumOfObs*5];

	double * vVol=new double [iNumOfObs];

//	fstream fs;
//	string sTest="XXX.csv";
//	int iNumOfTest=3;
//	int iMax=10;
//	double * vTest= new double [iNumOfTest];
//	for(int k=0; k<iNumOfTest;k++)
//	{
//		vTest[k]+=k;
//	}
//
//	for(int j=0; j<iMax; j++)
//	{
//		for(int k=0; k<iNumOfTest;k++)
//		{
//			vTest[k]+=j;
//		}
//
//		SequentialStore(j, iMax,iNumOfTest, vTest , fs,  sTest);
//		cout <<"seq"<<endl;
//	}
//
//	delete [] vTest;
	cout<<"start"<<endl;
	for( int i=0; i<iNumOfObs; i++)
		{
			cout<<"data "<<vData[i]<<endl;
		}
//	fstream Xfs;
//	fstream Sfs;
//	string sXDraw="vXDraw.csv";
//	string sSDraw="vSDraw.csv";
	for( int i=0; i<iNumOfIter; i++)
	{
		cout <<"==== Iteration: " << i << " ===="<< endl;
		/*Draw vY Star */
		DrawYStarT( sParam, vData,iNumOfObs ,  gsl_random_num);
		cout << "YStar is Done!" <<endl;
		/*Draw AuxY AuxH*/
		DrawAuxYandAuxHT(sParam, iNumOfObs , vData, gsl_random_num);
		cout<<"Aux Done "<<endl;

		/* Draw Gamma */
		DrawGammaAdaptiveRWT(vData,  iNumOfObs,i, gsl_random_num, sParam,dCovar,  dSum  );
		cout << "dGamma " << sParam.dGamma[0] << endl;

		/* Draw Nu */
		DrawNuAndLambdaT(sParam ,   iNumOfObs,  gsl_random_num,  0.1,  5);
		cout << "dNu " << sParam.dNu[0] << endl;
		/*Draw Phi, Sigma2*/
		DrawPhiSigmaAdaptiveRWWithSeasonal( iNumOfKnots, mWTilde,sParam, iNumOfObs,  i , gsl_random_num,  mCovar,   mSum );
		cout << "dPhi " << sParam.dPhi[0] << endl;
		cout << "dSigma " << sParam.dSigma2[0] << endl;

		/*Draw X, Mu*/

		DrawXandMuWithSeasonal(iNumOfKnots, mWTilde, sParam,iNumOfObs, gsl_random_num );
		cout << "dMu " << sParam.dMu[0] << endl;

		if(i>=iBurnIn)
		{
			for(int k=0; k<iNumOfObs; k++)
			{
				vVol[k]=sqrt((1-sParam.dGamma[0])*VarOrdT(exp((sParam.dMu[0]+sParam.vS[k]+sParam.vX[k])/2), sParam.dNu[0]));
			}

			if(i<iBurnIn+5)
			{
				for(int k=0; k<iNumOfObs; k++)
				{
					mMarkerXLow[k*5+i-iBurnIn]=sParam.vX[k];
					mMarkerSLow[k*5+i-iBurnIn]=sParam.vS[k];;
					mMarkerVolLow[k*5+i-iBurnIn]=vVol[k];

					mMarkerXHigh[k*5+i-iBurnIn]=sParam.vX[k];
					mMarkerSHigh[k*5+i-iBurnIn]=sParam.vS[k];;
					mMarkerVolHigh[k*5+i-iBurnIn]=vVol[k];
				}

			}
			else
			{
				if(i==iBurnIn+5)
				{
					InitializeQuantilePPAlgo(iNumOfObs, 0.05, mMarkerXLow,  mNXLow , mNDXLow);
					InitializeQuantilePPAlgo(iNumOfObs, 0.05, mMarkerSLow,  mNSLow , mNDSLow);
					InitializeQuantilePPAlgo(iNumOfObs, 0.05, mMarkerVolLow,  mNVolLow , mNDVolLow);
					InitializeQuantilePPAlgo(iNumOfObs, 0.95, mMarkerXHigh,  mNXHigh , mNDXHigh);
					InitializeQuantilePPAlgo(iNumOfObs, 0.95, mMarkerSHigh,  mNSHigh , mNDSHigh);
					InitializeQuantilePPAlgo(iNumOfObs, 0.95, mMarkerVolHigh,  mNVolHigh , mNDVolHigh);
				}

				 QuantilePPAlgo(iNumOfObs, 0.05,sParam.vX,  mMarkerXLow,  mNXLow , mNDXLow);
				 QuantilePPAlgo(iNumOfObs, 0.05, sParam.vS, mMarkerSLow,  mNSLow , mNDSLow);
				 QuantilePPAlgo(iNumOfObs, 0.05,vVol, mMarkerVolLow,  mNVolLow , mNDVolLow);
				 QuantilePPAlgo(iNumOfObs, 0.95,sParam.vX, mMarkerXHigh,  mNXHigh , mNDXHigh);
				 QuantilePPAlgo(iNumOfObs, 0.95,sParam.vS, mMarkerSHigh,  mNSHigh , mNDSHigh);
				 QuantilePPAlgo(iNumOfObs, 0.95,vVol, mMarkerVolHigh,  mNVolHigh , mNDVolHigh);
			}


//			SequentialStore(i-iBurnIn, iNumOfIter-iBurnIn,iNumOfObs, sParam.vX , Xfs,  sPrefix+sXDraw);
//			SequentialStore(i-iBurnIn, iNumOfIter-iBurnIn,iNumOfObs, sParam.vS , Sfs,  sPrefix+sSDraw);

			for(int j=0; j<iNumOfObs;j++)
			{
				vXEst[j]=vXEst[j]+(sParam.vX[j])/(iNumOfIter-iBurnIn);
				vSeasonEst[j]=vSeasonEst[j]+(sParam.vS[j])/(iNumOfIter-iBurnIn);
				vVolEst[j]=vVolEst[j]+vVol[j]/(iNumOfIter-iBurnIn);
				vLogIntEst[j]=vLogIntEst[j]+(sParam.dMu[0]+sParam.vS[j]+sParam.vX[j])/(iNumOfIter-iBurnIn);

				vYStarEst[j]=vYStarEst[j]+(sParam.vYStar[j])/(iNumOfIter-iBurnIn);
				vLambdaEst[j]=vLambdaEst[j]+(sParam.vLambda[j])/(iNumOfIter-iBurnIn);

				vX2Est[j]=vX2Est[j]+sParam.vX[j]*sParam.vX[j]/(iNumOfIter-iBurnIn);
				vSeason2Est[j]=vSeason2Est[j]+sParam.vS[j]*sParam.vS[j]/(iNumOfIter-iBurnIn);
				vVol2Est[j]=vVol2Est[j]+vVol[j]*vVol[j]/(iNumOfIter-iBurnIn);
			}


		}

		/*Save new draws */
		vector<double> vTempRow;

		vTempRow.push_back(sParam.dMu[0]);
		vTempRow.push_back(sParam.dPhi[0]);
		vTempRow.push_back(sParam.dSigma2[0]);
		vTempRow.push_back(sParam.dGamma[0]);
		for(int k=0; k<iNumOfKnots-1; k++)
		{
			vTempRow.push_back(sParam.vBeta[k]);
		}
		vTempRow.push_back(sParam.dNu[0]);


		mChain.push_back(vTempRow);
	}
	time (&end);
	double dif = difftime (end,start);
	printf ("Elasped time is %.2lf seconds. \n", dif );


		FILE * pFile;
		pFile = fopen ("OrdT_time.txt","w");
		fprintf (pFile, "OrdT time = %16.4f s\n",dif);
		fclose (pFile);
		
	for(int j=0; j<iNumOfObs;j++)
	{
		vXVolEst[j]=vX2Est[j]-vXEst[j]*vXEst[j];
		vSeasonVolEst[j]=vSeason2Est[j]-vSeasonEst[j]*vSeasonEst[j];
		vVolVolEst[j]=vVol2Est[j]-vVolEst[j]*vVolEst[j];
	}

//	string sXTrueFile="vXTrue.csv";
//	WriteOutDoubleArray(vXTrue, iNumOfObs, 1, sXTrueFile);

	string sXEstFile="vXEst.csv";
	WriteOutDoubleArray( vXEst , iNumOfObs, 1, sPrefix+sXEstFile);
	cout<<"X Est Saved!"<<endl;

//	string sLogIntTrueFile="vLogIntTrue.csv";
//	WriteOutDoubleArray(vLogIntTrue, iNumOfObs, 1, sLogIntTrueFile);

	string sLogIntEstFile="vLogIntEst.csv";
	WriteOutDoubleArray( vLogIntEst , iNumOfObs, 1, sPrefix+sLogIntEstFile);
	cout<<"LogInt Est Saved!"<<endl;

	string sSeasonEstFile="vSeasonEst.csv";
	WriteOutDoubleArray( vSeasonEst , iNumOfObs, 1, sPrefix+sSeasonEstFile);
	cout<<"Season Est Saved!"<<endl;


	string sYStarEstFile="vYStarEst.csv";
	WriteOutDoubleArray( vYStarEst , iNumOfObs, 1, sPrefix+sYStarEstFile);
	cout<<"YStar Est Saved!"<<endl;

//	string sYStarTrueFile="vYStarTrue.csv";
//	WriteOutDoubleArray( vYStarTrue , iNumOfObs, 1, sYStarTrueFile);
//
//	string sLambdaTrueFile="vLambdaTrue.csv";
//	WriteOutDoubleArray( vLambdaTrue , iNumOfObs, 1, sLambdaTrueFile);

	string sLambdaEstFile="vLambdaEst.csv";
	WriteOutDoubleArray(vLambdaEst , iNumOfObs, 1, sPrefix+sLambdaEstFile);
	cout<<"Lambda Est Saved!"<<endl;

	string sVolEstFile="vVolEst.csv";
	WriteOutDoubleArray( vVolEst , iNumOfObs, 1, sPrefix+sVolEstFile);
	cout<<"Vol Est Saved!"<<endl;

	string sXLowFile="vXLow.csv";
	WriteOutDoubleArray(  mMarkerXLow, iNumOfObs, 5, sPrefix+sXLowFile);
	cout<<"X Low Saved!"<<endl;

	string sXHighFile="vXHigh.csv";
	WriteOutDoubleArray(  mMarkerXHigh, iNumOfObs, 5, sPrefix+sXHighFile);
	cout<<"X High Saved!"<<endl;

	string sSLowFile="vSLow.csv";
	WriteOutDoubleArray(  mMarkerSLow, iNumOfObs, 5, sPrefix+sSLowFile);
	cout<<"S Low Saved!"<<endl;

	string sSHighFile="vSHigh.csv";
	WriteOutDoubleArray(  mMarkerSHigh, iNumOfObs, 5, sPrefix+sSHighFile);
	cout<<"S High Saved!"<<endl;

	string sVolLowFile="vVolLow.csv";
	WriteOutDoubleArray(  mMarkerVolLow, iNumOfObs, 5, sPrefix+sVolLowFile);
	cout<<"Vol Low Saved!"<<endl;

	string sVolHighFile="vVolHigh.csv";
	WriteOutDoubleArray(  mMarkerVolHigh, iNumOfObs, 5, sPrefix+sVolHighFile);
	cout<<"Vol High Saved!"<<endl;

	string sXVolFile="vXVol.csv";
	WriteOutDoubleArray(  vXVolEst, iNumOfObs, 1, sPrefix+sXVolFile);
	cout<<"X Vol Saved!"<<endl;

	string sSVolFile="vSVol.csv";
	WriteOutDoubleArray(  vSeasonVolEst, iNumOfObs, 1, sPrefix+sSVolFile);
	cout<<"S Vol Saved!"<<endl;

	string sVolVolFile="vVolVol.csv";
	WriteOutDoubleArray(  vVolVolEst, iNumOfObs, 1, sPrefix+sVolVolFile);
	cout<<"Vol Vol Saved!"<<endl;

	cout <<"Write out done"<<endl;


	delete [] sParam.vX;
	delete [] sParam.vBeta;
	delete [] sParam.vS;
	delete [] sParam.mAuxY;
	delete [] sParam.mAuxH;
	delete [] sParam.vYStar;
	delete [] sParam.vLambda;

	delete [] mMarkerXLow;
	delete [] mMarkerSLow;
	delete [] mMarkerVolLow;
	delete [] mMarkerXHigh;
	delete [] mMarkerSHigh;
	delete [] mMarkerVolHigh;

	delete [] mNXLow;
	delete [] mNSLow;
	delete [] mNVolLow;
	delete [] mNXHigh;
	delete [] mNSHigh;
	delete [] mNVolHigh;

	delete [] vVol;

	delete [] mNDXLow;
	delete [] mNDSLow;
	delete [] mNDVolLow;
	delete [] mNDXHigh;
	delete [] mNDSHigh;
	delete [] mNDVolHigh;

	delete [] vInitialLogVol;
	delete [] vLogIntEst;
//	delete [] vLogIntTrue;
//	delete [] mLogIntDraws;
//	delete [] mXDraws;
//	delete [] mNDraws;
//	delete [] mSeasonDraws;
//	delete [] mTau1Draws;
//	delete [] mTau2Draws;
	delete [] vXEst;
	delete [] vSeasonEst;
	delete [] vVolEst;
	delete [] vX2Est;
	delete [] vSeason2Est;
	delete [] vVol2Est;
	delete [] vYStarEst;
	delete [] vXVolEst;
	delete [] vSeasonVolEst;
	delete [] vVolVolEst;

	delete [] vKnots;

	delete dCovar;
	delete dSum;
	delete [] mCovar;
	delete [] mSum;

	delete sParam.iNumOfObs;
	delete sParam.dMu;
	delete sParam.dPhi;
	delete sParam.dSigma2;
	delete sParam.dGamma;
	delete sParam.dNu;
	delete sParam.dPriorNuA;
	delete sParam.dPriorNuB;
	delete sParam.dPriorGammaA;
	delete sParam.dPriorGammaB;
	delete sParam.dPriorMuMean;
	delete sParam.dPriorMuSigma2;
	delete sParam.dPriorPhiA;
	delete sParam.dPriorPhiB;
	delete sParam.dPriorSigmaA;
	delete sParam.dPriorSigmaB;

	cout <<"Clear memory done"<<endl;
}


//void MCMC::BootstrapFilterSkellam( vector<double> vParam, int iNumOfParticles)
//{
//	double dMu=vParam[0];
//	double dPhi=vParam[1];
//	double dSigma2=vParam[2];
//	double dGamma=vParam[3];
//
////	double dMaxLogWeight=-1e20;
//	double dSumWeights;
//	double dESS;
//
//	double* vParticles=new double[iNumOfParticles];
//	double* vLogWeights=new double[iNumOfParticles];
//	double* vWeights=new double[iNumOfParticles];
//	double* vNormalizedWeights=new double[iNumOfParticles];
//	double* vCumNormalizedWeights=new double[iNumOfParticles];
//
//	double* vMeans=new double[iNumOfObs];
//
//	const gsl_rng_type * T;
//	gsl_rng * gsl_random_num;
//	T=gsl_rng_default;
//	gsl_random_num=gsl_rng_alloc(T);
//	gsl_rng_set(gsl_random_num, 2);
//
//
//
//	for(int i=0; i<iNumOfObs; i++)
//	{
//		cout<<"XXXXXX Iter "<<i<<" data "<<vData[i] <<" XXXXXXXX"<<endl;
//		/* Initial step*/
//		if(i==0)
//		{
//			dSumWeights=0;
//			for(int j=0; j<iNumOfParticles; j++)
//			{
//				/* Initital draw */
//				vParticles[j]= DrawInital(dPhi ,dSigma2,  gsl_random_num);
////				cout << "init particles " << vParticles[j] <<endl;
//				/*Weights*/
//				vWeights[j]=ZeroSkellamPdf(vData[i],  dGamma, exp(dMu+vParticles[j]));
////				cout << "weights " << vWeights[j]<<endl;
//				/* Sum of weights */
//				dSumWeights+=vWeights[j];
//
//			}
//
//			/* Calculate statistics */
//			dESS=0;
//			vMeans[i]=0;
//			/* Normalized weights */
//			for(int j=0; j<iNumOfParticles; j++)
//			{
//				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
//				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
////				cout << "normalized weights " << vNormalizedWeights[j]<<endl;
//				if(j==0)
//				{
//					vCumNormalizedWeights[j]=vNormalizedWeights[j];
//				}
//				else
//				{
//					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
//				}
////				cout << "cum normalized weights " << vCumNormalizedWeights[j]<<endl;
//				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
//			}
//			dESS=1/dESS;
//
////			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
//		}
//		else
//		{
////			cout << "XXX ESS "<<dESS<<" XXX"<<endl;
//			/*Resample */
////			if(dESS<0.5*iNumOfParticles)
////			{
//
////				cout << "XXX Resampling  XXX"<<endl;
////				for(int k=0; k<iNumOfParticles; k++)
////				{
////					cout<<"Old Partciles " << vParticles[k] <<" Weights " << vNormalizedWeights[k] <<endl;
////				}
//				SystemicResampling(vParticles, vCumNormalizedWeights,iNumOfParticles,gsl_random_num);
////				for(int k=0; k<iNumOfParticles; k++)
////				{
////					cout<<"New Partciles " << vParticles[k] <<endl;
////				}
////			}
//
//			dSumWeights=0;
//			for(int j=0; j<iNumOfParticles; j++)
//			{
//				/* Propagate */
//				vParticles[j]=DrawTransition(vParticles[j],dPhi ,dSigma2,  gsl_random_num);
////				cout << "particles " << vParticles[j] <<endl;
//				/*Weights*/
////				if(dESS<0.5*iNumOfParticles)
////				{
//					vWeights[j]=ZeroSkellamPdf(vData[i],  dGamma, exp(dMu+vParticles[j]))/iNumOfParticles;
////				}
////				else
////				{
////					vWeights[j]=ZeroSkellamPdf(vData[i],  dGamma, exp(dMu+vParticles[j]))*vWeights[j];
////				}
////				cout << "weights " << vWeights[j]<<endl;
//				/* Sum of weights */
//				dSumWeights+=vWeights[j];
//			}
//
//
//
//
//			dESS=0;
//			vMeans[i]=0;
//			for(int j=0; j<iNumOfParticles; j++)
//			{
//
//				/* Normalized weights */
//				vNormalizedWeights[j]=vWeights[j]/dSumWeights;
//				/* Calculate statistics */
//				vMeans[i]+=vNormalizedWeights[j]*vParticles[j];
////				cout << "normalized weights at " << j <<" is "<< vNormalizedWeights[j]<<endl;
//				if(j==0)
//				{
////					cout<< "first "<<vNormalizedWeights[0] <<endl;
//					vCumNormalizedWeights[j]=vNormalizedWeights[j];
//				}
//				else
//				{
////					cout << "prev cum normalized weights at " << j-1 <<" is "<<  vCumNormalizedWeights[j-1]<<endl;
//					vCumNormalizedWeights[j]=vCumNormalizedWeights[j-1]+vNormalizedWeights[j];
//				}
////				cout << "cum normalized weights at " << j <<" is "<<  vCumNormalizedWeights[j]<<endl;
//				dESS+=vNormalizedWeights[j]*vNormalizedWeights[j];
//			}
//			dESS=1/dESS;
//		}
//
//
//	}
//
//	string sXTrue="vXTrue.csv";
//	WriteOutDoubleArray( vXTrue , iNumOfObs, 1, sXTrue);
//
//	string sBFMeans="vXMeansBF.csv";
//	WriteOutDoubleArray( vMeans , iNumOfObs, 1, sBFMeans);
//
//	delete [] vParticles;
//	delete [] vMeans;
//	delete [] vLogWeights;
//	delete [] vNormalizedWeights;
//	delete [] vWeights;
//}
