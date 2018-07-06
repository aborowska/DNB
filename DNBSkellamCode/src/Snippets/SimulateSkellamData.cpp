
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
void MCMC::SimulateSkellamData( int iNumOfSim)
{
	string sType;
	sType="Sk";
	const gsl_rng_type * T;
	gsl_rng * r;

	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, 1);

	vData = new double[iNumOfSim];
	vXTrue = new double[iNumOfSim];
	vNTrue=new double[iNumOfSim];
	vTau1True=new double[iNumOfSim];
	vTau2True=new double[iNumOfSim];

	vSeasonTrue = new double[iNumOfSim];
	vTimes = new double[iNumOfSim];
	iNumOfObs=iNumOfSim;

	cout<<"dMuTrue "<<dMuTrue<<endl;
	cout<<"dPhiTrue "<<dPhiTrue<<endl;
	cout<<"dSigma2True "<<dSigma2True<<endl;
	cout<<"dGammaTrue "<<dGammaTrue<<endl;

	double dStateSigma=sqrt(dSigma2True);

	int iNumOfSimPerDay=2000;
	// SimulateSeasonal(iNumOfSimPerDay, iNumOfSim/iNumOfSimPerDay,  vTimes,  vSeasonTrue);
	SimulateSeasonal_str(sType, iNumOfSimPerDay, iNumOfSim/iNumOfSimPerDay,  vTimes,  vSeasonTrue);

	for(int i=0; i<iNumOfSim; i++){
		if(i==0)
		{
			/* Log intensity */
			vXTrue[i]= gsl_ran_gaussian(r,sqrt(dSigma2True/(1-dPhiTrue*dPhiTrue)) );
		}
		else
		{
			vXTrue [i]= dPhiTrue * vXTrue[i-1]+ gsl_ran_gaussian(r,dStateSigma) ;
		}
//		vXTrue [i]=0;
		/* Conditional Skellam return */
		double dPoissionMean= exp(dMuTrue +vSeasonTrue[i] + vXTrue[i]);
//		int Poission1=gsl_ran_poisson(r,dPoissionMean);
//		int Poission2= gsl_ran_poisson(r,dPoissionMean);
//		int iSkellam=Poission1 -  Poission2;
//		vNTrue[i]=Poission1 +  Poission2;

		int iSkellam=0;
		int iN=0;
		double dTau1=0;
		double dTau2=0;
		SimulateSkellam(dPoissionMean, r, &iSkellam, &iN, &dTau1, &dTau2);

		vNTrue[i]=iN;
		vTau1True[i]=dTau1;
		vTau2True[i]=dTau2;

		double dU=gsl_rng_uniform(r);
		if(dU<=dGammaTrue)
		{
			vData[i]=0;
		}
		else
		{
			vData[i]= iSkellam ;
		}
		cout << "Data at "<< i<<" is " << vData[i]<<  endl;
	}
	gsl_rng_free(r);
	return;
}



void SimulateSkellam(double dLambda, const gsl_rng * gsl_random_num, int*  iY, int* iN, double* dTau1, double * dTau2)
{
	double dTau=0;
	double dTauPrev=0;

	int iNTemp=0;
	int iYTemp=0;

	int iNPrevTemp=0;
	int iYPrevTemp=0;

	while(dTau < 1)
	{
		iNPrevTemp=iNTemp;
		iYPrevTemp=iYTemp;
		dTauPrev=dTau;
		dTau=dTau+gsl_ran_exponential(gsl_random_num, 1/(2*dLambda));

		iNTemp=iNTemp+1;
		double dU=gsl_rng_uniform(gsl_random_num);
		if(dU<=0.5)
		{
			iYTemp=iYTemp+1;
		}
		else
		{
			iYTemp=iYTemp-1;
		}
	}

	iN[0]=iNPrevTemp;
	iY[0]=iYPrevTemp;
	dTau1[0]=dTau-dTauPrev;
	dTau2[0]=dTauPrev;
	cout <<"while end"<< endl;
}

