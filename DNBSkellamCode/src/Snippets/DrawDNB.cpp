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


void DrawNuInBlockDiscreteDNB(double * vData, struct AllParam sParam ,int iNumOfObs, const gsl_rng * gsl_random_num, double dResolution,  double dDelta)
		// DrawNuInBlockDiscreteDNB(vData, sParam ,  iNumOfObs, gsl_random_num, 0.2,3);

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





void NuInBlockLogPosteriorDNB(double * dLogNu, double *vData, struct AllParam sParam,  double *dPriorA, double *dPriorB, unsigned int iNumOfObs, double * dLogPosterior)
{

	double dNu=exp(dLogNu[0]);
	double dOldNu=sParam.dNu[0];
	sParam.dNu[0]=dNu;

	CaculateDNBLL(vData, iNumOfObs, sParam,dLogPosterior );

	dLogPosterior[0]+=dPriorA[0]*log(dPriorB[0]) +dPriorA[0]*log(dNu) -dPriorB[0]*dNu-lgamma(dPriorA[0]);
	sParam.dNu[0]=dOldNu;
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


/* ********************************************* */
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

/* ********************************************* */
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

/* ********************************************* */
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

//		cout << "CumU " << dCumU[0] << " dU "<< dU[0]<< endl;
		while(dCumU[0]<dU[0])
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

/* ********************************************* */
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

/* ********************************************* */
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
			/* Draw indicator form the mixture approximation for tau2*/
//			for(int i=0 ; i<10; i++)
//			{
//				cout<< "vWeight at "<<i<<" is " << vWeights[i]<<endl;
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


/* ********************************************* */
void DrawGammaAdaptiveRWDNB(  int iNumOfObs,  int iNumOfIter, const gsl_rng * gsl_random_num, struct AllParam sParam,
		double* vIndicator, double* vDNB, double* dCovar, double* dSum, double* dM2 )
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
/* 	else if(iNumOfIter==2)
	{
		dSigma=pow(dCovar[0]-0.5*dSum[0]*dSum[0],0.5);
		if(dSigma==0)
		{
			dSigma=0.1;
		}
	} */
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
//	if(iNumOfIter<=1)
// >>>>>>>>>>>>>> dSum is the MEAN now	<<<<<<<<<<<<<<<<<<
	if(iNumOfIter<1)
	{
/* 		dCovar[0]=dCovar[0]+pow(LogitTransform(sParam.dGamma[0]),2);
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]); */
		dSum[0] = dSum[0]+LogitTransform(sParam.dGamma[0]);
		dM2[0] = 0;
		dCovar[0] = 0;		
	}
	// else if(iNumOfIter==2)
	// {
// /* 		dCovar[0]=(iNumOfIter-1)*pow(dSigma,2)/iNumOfIter+dSum[0]*dSum[0]/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		// dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		// dCovar[0]=dCovar[0]-dSum[0]*dSum[0]/(iNumOfIter*(iNumOfIter+1)); */
		// double dSum_old = dSum[0];
		// dSum[0] = dSum[0]+LogitTransform(sParam.dGamma[0]);
		// dM2[0] = dM2[0] + (LogitTransform(sParam.dGamma[0]) - dSum_old/iNumOfIter)*(LogitTransform(sParam.dGamma[0]) - dSum[0]/(iNumOfIter+1));		
		// dCovar[0] = dM2[0]/iNumOfIter;		
	// }
	else
	{
/* 		dCovar[0]=(iNumOfIter-1)*dCovar[0]/iNumOfIter+(dSum[0]/iNumOfIter)*(dSum[0]/iNumOfIter)+pow(LogitTransform(sParam.dGamma[0]),2)/iNumOfIter;
		dSum[0]=dSum[0]+LogitTransform(sParam.dGamma[0]);
		dCovar[0]=dCovar[0]-(dSum[0]/iNumOfIter)*(dSum[0]/(iNumOfIter+1)); */	
		double dSum_old = dSum[0];
		dSum[0] = dSum[0] + (LogitTransform(sParam.dGamma[0]) - dSum_old)/(iNumOfIter+1);
		dM2[0] = dM2[0] + (LogitTransform(sParam.dGamma[0]) - dSum_old)*(LogitTransform(sParam.dGamma[0]) - dSum[0]);
		dCovar[0] = dM2[0]/iNumOfIter;
		
	}

	cout << ":: ARW3 dMean " << dSum[0]<< endl;
	cout << ":: ARW3 dSigma " << dSigma<< endl;
	cout << ":: ARW3 dCovar " << dCovar[0]<< endl;
	cout << ":: ARW proposed logit " << dNewGammaLogit << "proposed  " << LogitTransformBack(dNewGammaLogit)<< endl;
	cout << ":: ARW accepted logit " << LogitTransform(sParam.dGamma[0])<< "accepted " << sParam.dGamma[0] <<  endl;

	delete dLogPosteriorNew;
	delete dLogPosteriorOld;

}

/* ********************************************* */
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

/* ********************************************* */


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

/* ********************************************* */


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

/* ********************************************* */

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


/* ********************************************* */

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
 
	/* AHat */
	KalmanFilterWithSeasonal(iNumOfKnots,   mWTilde ,mY, sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL);
	KalmanSmootherWithSeasonal( iNumOfKnots,    mWTilde ,sParam.vN,mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHat,  mVHat );

	/* AHatPlus */

	KalmanFilterWithSeasonal(iNumOfKnots,   mWTilde ,mYPlus, sParam.vN, mH, mZ,  mT, mQ, mA1,  mP1, iDimOfObs,iDimOfStates,  iNumOfObs, mA,  mP,  mFInv,  mV,  mL); 
	KalmanSmootherWithSeasonal( iNumOfKnots,    mWTilde ,sParam.vN, mZ,  mA, mP, mV,  mFInv,  mL,   iDimOfObs, iDimOfStates, iNumOfObs, mAHatPlus,  mVHatPlus);
 

	double * mAHatAHatPlus =new double[iDimOfStates*iNumOfObs];
	MatrixSub(mAHat, mAHatPlus, iNumOfObs, iDimOfStates, mAHatAHatPlus);
	MatrixAdd(mAHatAHatPlus, mAPlus, iNumOfObs, iDimOfStates, mDraw); 

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

/* ********************************************* */


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

/* ********************************************* */


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



/* ********************************************* */




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
