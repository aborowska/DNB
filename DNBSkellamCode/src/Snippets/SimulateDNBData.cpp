/*
 * Simulate  from dynamic zero inflated negative binomial difference model
 *
 * iNumOfSim	number of simulated data points
 *
 * Note that the algorithm uses :
 * dMuTrue;
 * dPhiTrue;
 * Sigma2True;
 * dGammaTrue;
 * dNuTrue;
 *
 */
void MCMC::SimulateDNBData( int iNumOfSim)
{
	string sType;
	sType="DNB";
	const gsl_rng_type * T;
	gsl_rng * r;

	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, iSeedNo);

	vData = new double[iNumOfSim];
	vXTrue  = new double[iNumOfSim];

	vSeasonTrue = new double[iNumOfSim];
	vTimes = new double[iNumOfSim];
	iNumOfObs=iNumOfSim;

	cout<<"dMuTrue "<< dMuTrue<<endl;
	cout<<"dPhiTrue "<< dPhiTrue<<endl;
	cout<<"dSigma2True "<< dSigma2True <<endl;
	cout<<"dGammaTrue "<< dGammaTrue <<endl;
	cout<<"dNuTrue "<< dNuTrue<<endl;
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
		/* Conditional DNB return */
		double dPoissionMean1= exp(dMuTrue +vSeasonTrue[i]+ vXTrue[i]) * gsl_ran_gamma(r, dNuTrue, 1/dNuTrue);
		double dPoissionMean2= exp(dMuTrue +vSeasonTrue[i]+vXTrue[i]) * gsl_ran_gamma(r, dNuTrue, 1/dNuTrue);
		int iDNB =gsl_ran_poisson(r,dPoissionMean1) -  gsl_ran_poisson(r,dPoissionMean2) ;
		double dU=gsl_rng_uniform(r);
		if(dU<=dGammaTrue)
		{
			vData[i]=0;
		}
		else
		{
			vData[i]=iDNB;
		}

		cout << "Data at "<< i<<" is " << vData[i]<<  endl;
	}

	gsl_rng_free(r);
	return;
}


void SimulateSeasonal_str(string sType, int iNumOfSimPerDay, int iNumOfDays, double* vTimes, double * vSeason)
{
	string sPrefix;
	string sInput;
	CreatPrefix(sInput,sType, &sPrefix );
		
	int iTimeInSec=23400;
	double dDuration=(double)iTimeInSec/iNumOfSimPerDay;

	double *vUnnormalizedSeason=new double[iNumOfSimPerDay*iNumOfDays];
	double *dDailyAvgSeason=new double;

	/* Unnormalized sesonal pattern */
	for(int d=0; d<iNumOfDays; d++)
	{
		dDailyAvgSeason[0]=0;
		for(int i=0; i<iNumOfSimPerDay; i++)
		{
			vTimes[i+d*iNumOfSimPerDay]=(i+1)*dDuration;
			vUnnormalizedSeason[i+d*iNumOfSimPerDay]=(vTimes[i]-10800)*(vTimes[i]-10800)/(58320000);
			dDailyAvgSeason[0]=dDailyAvgSeason[0]+vUnnormalizedSeason[i];
		}

		dDailyAvgSeason[0]=dDailyAvgSeason[0]/iNumOfSimPerDay;

		/* Zero mean normalized seasonal pattern */
		for(int i=0; i<iNumOfSimPerDay; i++)
		{
			vSeason[i+d*iNumOfSimPerDay]=vUnnormalizedSeason[i+d*iNumOfSimPerDay]-dDailyAvgSeason[0];
//			vSeason[i+d*iNumOfSimPerDay]=0;
		}
	}

	/* Saving results */
	string sSeasonFile=(sPrefix+"vSeason.csv").c_str();
	WriteOutDoubleArray(vSeason, iNumOfSimPerDay*iNumOfDays, 1, sSeasonFile);

	string sTimeFile=(sPrefix+"vTime.csv").c_str();
	WriteOutDoubleArray(vTimes, iNumOfSimPerDay*iNumOfDays, 1, sTimeFile);

	delete [] vUnnormalizedSeason;
	delete dDailyAvgSeason;
}
