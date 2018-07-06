void MCMC::EstimateSkellam( int iNumOfIter, int iBurnIn)
{
	sType="Sk";
	CreatPrefix(sInput,sType, &sPrefix );
	cout<< "prefix done" <<endl;
	/* Initialize parameters and latent variables */
	const gsl_rng_type * T;
	gsl_rng * gsl_random_num;
	T=gsl_rng_default;
	gsl_random_num=gsl_rng_alloc(T);
	gsl_rng_set(gsl_random_num, 41);

	double* vZeroSkellamDens = new double[iNumOfObs];
	double* vSkellamDens = new double[iNumOfObs];
	double* vLogBessel = new double[iNumOfObs];
	double* vIndicator = new double[iNumOfObs];

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
	sParam.vN = new  int  [iNumOfObs];
	sParam.vTau1= new double[iNumOfObs];
	sParam.vTau2= new double[iNumOfObs];
	sParam.mAuxY=new double[2*iNumOfObs];
	sParam.mAuxH=new double[4*iNumOfObs];

	sParam.iNumOfObs[0]=iNumOfObs;

	sParam.dMu[0]=0;
	sParam.dPhi[0]=0.95;
	sParam.dSigma2[0]=0.01;
	sParam.dGamma[0]=0.3;
//	sParam.dGamma[0]=0;

//	sParam.dMu[0]=2;
//	sParam.dPhi[0]=0.90;
//	sParam.dSigma2[0]=0.5;
//	sParam.dGamma[0]=0.2;

	sParam.dPriorGammaA[0]=1.7;
	sParam.dPriorGammaB[0]=10;
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
//		sParam.vS[i]=vSeasonTrue[i];
//		sParam.vX[i]=vXTrue[i];
//		sParam.vS[i]=vSeasonTrue[i];
//		sParam.vX[i]=vXTrue[i];
//		sParam.vN[i]=vNTrue[i];
//		sParam.vTau1[i]=vTau1True[i];
//		sParam.vTau2[i]=vTau2True[i];
	}

	/*Initializing the iteration*/
	time_t start,end;
	time (&start);

	double * vLogIntEst=new double[iNumOfObs];
	double * vXEst=new double[iNumOfObs];
	double * vSeasonEst=new double[iNumOfObs];
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
		vLambdaEst[j]=0;
		vXVolEst[j]=0;
		vSeasonVolEst[j]=0;
		vVolVolEst[j]=0;
	}

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

//	string svInitLogIntFile="vInitialLogVol.csv";
//	WriteOutDoubleArray( vInitialLogVol, iNumOfObs, 1, svInitLogIntFile);
//	string svInitSFile="vInitS.csv";
//	WriteOutDoubleArray( sParam.vS, iNumOfObs, 1, svInitSFile);
//	string svInitXFile="vInitX.csv";
//	WriteOutDoubleArray( sParam.vX, iNumOfObs, 1, svInitXFile);

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
		CaculateSkellam(vData,iNumOfObs, sParam, vZeroSkellamDens,vIndicator, vSkellamDens,vLogBessel);
		cout<<"Skellam Done "<<endl;
//		/*Draw N, Tau1, Tau2, S */
//		DrawAuxVariables(i, vData, iNumOfObs, sParam, vZeroSkellamDens);
		/*Draw N */
		DrawN(  sParam,  vData,  iNumOfObs ,  gsl_random_num, vZeroSkellamDens,vLogBessel );
		cout<<"N Done "<<endl;
		/*Draw Tau1, Tau2 */
		DrawTau( sParam, iNumOfObs ,  gsl_random_num);
		cout<<"Tau Done "<<endl;
		/*Draw AuxY AuxH*/
		DrawAuxYandAuxH(sParam, iNumOfObs , gsl_random_num);
		cout<<"Aux Done "<<endl;

		/* Draw Gamma */
//		DrawGammaLaplaceApprox(iNumOfObs, gsl_random_num, sParam,vIndicator, vSkellamDens  );
		DrawGammaAdaptiveRW(  iNumOfObs,i, gsl_random_num, sParam,vIndicator, vSkellamDens, dCovar,  dSum  );

		cout << "dGamma " << sParam.dGamma[0] << endl;

		/*Draw Seasonal*/
//		DrawSesonal(sParam, iNumOfObs, iNumOfKnots,mWTilde ,gsl_random_num);

		/*Draw Phi, Sigma2*/
//		DrawPhiSigmaFullConditional( sParam, iNumOfObs,   gsl_random_num );
//		DrawPhiSigmaLaplaceApprox(  sParam, iNumOfObs,   gsl_random_num );
//		DrawPhiSigmaLaplaceApproxWithSeasonal(  sParam, iNumOfObs,   gsl_random_num );
		DrawPhiSigmaAdaptiveRWWithSeasonal( iNumOfKnots, mWTilde,sParam, iNumOfObs,  i , gsl_random_num,  mCovar,   mSum );
//		DrawPhiSigmaLaplaceApproxSGD(  sParam, iNumOfObs,   gsl_random_num );
//		DrawPhiSigmaAdaptiveRW( sParam, iNumOfObs,  i , gsl_random_num,  mCovar,   mSum );
		cout << "dPhi " << sParam.dPhi[0] << endl;
		cout << "dSigma " << sParam.dSigma2[0] << endl;

		/*Draw X, Mu*/
//		DrawXandMu(sParam,iNumOfObs, gsl_random_num );
		DrawXandMuWithSeasonal(iNumOfKnots, mWTilde, sParam,iNumOfObs, gsl_random_num );
		cout << "dMu " << sParam.dMu[0] << endl;

		if(i>=iBurnIn)
		{
			for(int k=0; k<iNumOfObs; k++)
			{
				vVol[k]=sqrt((1-sParam.dGamma[0])*2*exp((sParam.dMu[0]+sParam.vS[k]+sParam.vX[k])));
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
	printf ("Elapsed time is %.2lf seconds. \n", dif );

		FILE * pFile;
		pFile = fopen ("Sk_time.txt","w");
		fprintf (pFile, "Skellman time = %16.4f s\n",dif);
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

//	string sYStarTrueFile="vYStarTrue.csv";
//	WriteOutDoubleArray( vYStarTrue , iNumOfObs, 1, sYStarTrueFile);
//
//	string sLambdaTrueFile="vLambdaTrue.csv";
//	WriteOutDoubleArray( vLambdaTrue , iNumOfObs, 1, sLambdaTrueFile);


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

	cout <<"Write out done"<<endl;

	delete [] vZeroSkellamDens;
	delete [] vIndicator;
	delete [] vLogBessel;
	delete [] vSkellamDens;
	delete [] sParam.vX;
	delete [] sParam.vBeta;
	delete [] sParam.vS;
	delete [] sParam.vN;
	delete [] sParam.vTau1;
	delete [] sParam.vTau2;
	delete [] sParam.mAuxY;
	delete [] sParam.mAuxH;


	delete [] vLogIntEst;
	delete [] vXEst;
	delete [] vSeasonEst;
	delete [] vVolEst;
	delete [] vX2Est;
	delete [] vSeason2Est;
	delete [] vVol2Est;
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
