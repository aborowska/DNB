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

//	cout << "dSum " << dSum[0]<< endl;
//	cout << "dCovar " << dCovar[0]<< endl;
//	cout << "proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
//	cout << "accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;


}
