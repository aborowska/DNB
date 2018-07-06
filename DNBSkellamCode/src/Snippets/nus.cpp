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



void NuInBlockLogPosteriorDNB(double * dLogNu, double *vData, struct AllParam sParam,  double *dPriorA, double *dPriorB, unsigned int iNumOfObs,
	double * dLogPosterior)
{

	double dNu=exp(dLogNu[0]);
	double dOldNu=sParam.dNu[0];
	sParam.dNu[0]=dNu;



	CaculateDNBLL(vData, iNumOfObs, sParam,dLogPosterior );


	dLogPosterior[0]+=dPriorA[0]*log(dPriorB[0]) +dPriorA[0]*log(dNu) -dPriorB[0]*dNu-lgamma(dPriorA[0]);
	sParam.dNu[0]=dOldNu;

}


void DrawNuInBlockDiscreteDNB(double * vData, struct AllParam sParam ,int iNumOfObs, const gsl_rng * gsl_random_num, double dResolution,  double dDelta)
// double dResolution 0.2,  double dDelta 3
{
	/* Propose draw */
	double dU=gsl_rng_uniform(gsl_random_num);
	int iN=2*dDelta/dResolution+1; // = 31
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

	if((dNuProposed <= 128) & (dNuProposed >=0.1))
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
