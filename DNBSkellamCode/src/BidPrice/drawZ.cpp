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
