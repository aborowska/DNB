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


void CaculateDNB(const double * vData, int iNumOfObs,  struct AllParam sParam,double * vDNB, double * vZeroDNBDens, double * vIndicator)
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
