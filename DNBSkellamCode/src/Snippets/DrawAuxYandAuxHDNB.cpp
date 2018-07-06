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
//   		cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
//			cout<< "vVar at "<<i<<" is " << vVariances[i]<<endl;
//		}
//

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
///				cout<< "vMean at "<<i<<" is " << vMeans[i]<<endl;
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