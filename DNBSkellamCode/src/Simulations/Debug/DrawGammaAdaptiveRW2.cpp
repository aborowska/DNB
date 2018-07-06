void DrawGammaAdaptiveRW(  int iNumOfObs, int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam, double* vIndicator, double* vSkellam, double* dCovar,  double* dSum  )
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
    
    double dAccept=fmin(dLogPosteriorNew[0]-dLogPosteriorOld[0],0);
	// double dAccept=fmin(-GammaLogPosterior(newx, &sObjectiveFunParam )*iNumOfObs+GammaLogPosterior(oldx, &sObjectiveFunParam)*iNumOfObs,0);

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

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;
}
