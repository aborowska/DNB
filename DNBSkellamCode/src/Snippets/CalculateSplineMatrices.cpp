void CalculateSplineMatrices( int iNumOfObs, double* vTimes,  int iNumOfKnots, double* vKnots, gsl_matrix * mWTilde)
{
	gsl_matrix * mTheta = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_matrix * mLambda = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_matrix * mP = gsl_matrix_alloc (iNumOfObs, iNumOfKnots);
	gsl_matrix * mW = gsl_matrix_alloc (iNumOfObs, iNumOfKnots);
	gsl_matrix * mWStar = gsl_matrix_alloc (1, iNumOfKnots);


	double *dLambda=new double;
	double *dH=new double;
	double *dHPlus=new double;

	/* Lambda and Theta */
	for(int i=0; i<iNumOfKnots; i++)
	{
		for(int j=0; j<iNumOfKnots; j++)
		{
			gsl_matrix_set (mLambda, i, j, 0);
			gsl_matrix_set (mTheta, i, j, 0);
		}

		if(i==0 ||  i==iNumOfKnots-1 )
		{
			gsl_matrix_set (mLambda , i, i, 2);
		}
		else
		{
			dH[0]=vKnots[i]-vKnots[i-1];
			dHPlus[0]=vKnots[i+1]-vKnots[i];
			dLambda[0]=dHPlus[0]/(dH[0]+dHPlus[0]);
			/* Lambda matrix */
			gsl_matrix_set (mLambda , i, i-1, 1-dLambda[0]);
			gsl_matrix_set (mLambda, i, i, 2);
			gsl_matrix_set (mLambda , i, i+1, dLambda[0]);
			/* Theta matrix */
			gsl_matrix_set (mTheta , i, i-1, 6/(dH[0]*(dH[0]+dHPlus[0])));
			gsl_matrix_set (mTheta , i, i, -6/(dH[0]*dHPlus[0]));
			gsl_matrix_set (mTheta , i, i+1,6/(dHPlus[0]*(dH[0]+dHPlus[0])));
		}
	}

	/* P and Q */
	for(int i=0; i<iNumOfObs;i++)
	{
		int iTempKnot=1;
		while(vTimes[i] > vKnots[iTempKnot])
		{
//			cout<<"vTimes "<<vTimes[i]<<"vKnots "<<vKnots[iTempKnot]<<endl;
			iTempKnot=iTempKnot+1;
		}
//		cout <<"END "<<endl;
		for(int j=0; j<iNumOfKnots; j++)
		{
			if(j==iTempKnot-1)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mP , i, j, (vKnots[iTempKnot]-vTimes[i])*(pow(vKnots[iTempKnot]-vTimes[i],2)-dH[0]*dH[0])/(6*dH[0]) );
				gsl_matrix_set (mW , i, j, (vKnots[iTempKnot]-vTimes[i])/dH[0] );

			}
			else if (j==iTempKnot)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mP , i, j, (vTimes[i]-vKnots[iTempKnot-1])*(pow(vTimes[i]-vKnots[iTempKnot-1],2)-dH[0]*dH[0])/(6*dH[0]));
				gsl_matrix_set (mW , i, j, (vTimes[i]-vKnots[iTempKnot-1])/dH[0]);
			}
			else
			{
				gsl_matrix_set (mP , i, j, 0);
				gsl_matrix_set (mW , i, j, 0);
			}
		}
	}

	/* Lambda  inverse */
	gsl_matrix * mLambdaInv = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots);

	int s;
	gsl_linalg_LU_decomp (mLambda, mPermutation, &s);
	gsl_linalg_LU_invert (mLambda, mPermutation, mLambdaInv);


	/* W matrix */
	gsl_matrix * mLambdaInvTheta = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	for(int i=0; i<iNumOfKnots; i++)
	{
		for(int j=0; j<iNumOfKnots; j++)
		{
			gsl_matrix_set ( mLambdaInvTheta  ,1,j,0 );
		}
	}

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mLambdaInv,mTheta,0.0, mLambdaInvTheta);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mP,mLambdaInvTheta,1.0, mW);


	/* WStar matrix */
//	cout <<"calculate star started "<<endl;
	CalculateWStar( iNumOfObs, vTimes, iNumOfKnots,  vKnots, mLambdaInvTheta, mWStar );
//	cout <<"calculate star finished "<<endl;


	/* WTilde matrix */
	for(int i=0; i<iNumOfObs; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			gsl_matrix_set (mWTilde ,i,j, gsl_matrix_get (mW ,i,j) - gsl_matrix_get (mW ,i,iNumOfKnots-1)*gsl_matrix_get (mWStar ,0,j)/gsl_matrix_get (mWStar ,0,iNumOfKnots-1) );
		}
	}


	gsl_matrix_free (mTheta);
	gsl_matrix_free (mLambda);
	gsl_matrix_free (mP);
	gsl_matrix_free (mW);
	gsl_matrix_free (mWStar);
	gsl_matrix_free (mLambdaInv);

	gsl_permutation_free (mPermutation);

	delete dLambda;
	delete dH;
	delete dHPlus;

}





void CalculateSplineMatricesOut(int iNumOfObs, double* vTimes,  int iNumOfKnots, double* vKnots, gsl_matrix * mWTilde)
{
	gsl_matrix * mTheta = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_matrix * mLambda = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_matrix * mP = gsl_matrix_alloc (iNumOfObs, iNumOfKnots);
	gsl_matrix * mW = gsl_matrix_alloc (iNumOfObs, iNumOfKnots);
	gsl_matrix * mWStar = gsl_matrix_alloc (1, iNumOfKnots);

	double *dLambda=new double;
	double *dH=new double;
	double *dHPlus=new double;
	/* Lambda and Theta */
	for(int i=0; i<iNumOfKnots; i++)
	{
		for(int j=0; j<iNumOfKnots; j++)
		{
			gsl_matrix_set (mLambda, i, j, 0);
			gsl_matrix_set (mTheta, i, j, 0);
		}

		if(i==0 ||  i==iNumOfKnots-1 )
		{
			gsl_matrix_set (mLambda , i, i, 2);
		}
		else
		{
			dH[0]=vKnots[i]-vKnots[i-1];
			dHPlus[0]=vKnots[i+1]-vKnots[i];
			dLambda[0]=dHPlus[0]/(dH[0]+dHPlus[0]);
			/* Lambda matrix */
			gsl_matrix_set (mLambda , i, i-1, 1-dLambda[0]);
			gsl_matrix_set (mLambda, i, i, 2);
			gsl_matrix_set (mLambda , i, i+1, dLambda[0]);
			/* Theta matrix */
			gsl_matrix_set (mTheta , i, i-1, 6/(dH[0]*(dH[0]+dHPlus[0])));
			gsl_matrix_set (mTheta , i, i, -6/(dH[0]*dHPlus[0]));
			gsl_matrix_set (mTheta , i, i+1,6/(dHPlus[0]*(dH[0]+dHPlus[0])));
		}

	}

	/* P and Q */
	for(int i=0; i<iNumOfObs;i++)
	{
		int iTempKnot=1;
		while(vTimes[i] > vKnots[iTempKnot])
		{
//			cout<<"vTimes "<<vTimes[i]<<"vKnots "<<vKnots[iTempKnot]<<endl;
			iTempKnot=iTempKnot+1;
		}
//		cout <<"END "<<endl;
		for(int j=0; j<iNumOfKnots; j++)
		{
			if(j==iTempKnot-1)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mP , i, j, (vKnots[iTempKnot]-vTimes[i])*(pow(vKnots[iTempKnot]-vTimes[i],2)-dH[0]*dH[0])/(6*dH[0]) );
				gsl_matrix_set (mW , i, j, (vKnots[iTempKnot]-vTimes[i])/dH[0] );
			}
			else if (j==iTempKnot)
			{
				dH[0]=vKnots[iTempKnot]-vKnots[iTempKnot-1];
				gsl_matrix_set (mP , i, j, (vTimes[i]-vKnots[iTempKnot-1])*(pow(vTimes[i]-vKnots[iTempKnot-1],2)-dH[0]*dH[0])/(6*dH[0]));
				gsl_matrix_set (mW , i, j, (vTimes[i]-vKnots[iTempKnot-1])/dH[0]);
			}
			else
			{
				gsl_matrix_set (mP , i, j, 0);
				gsl_matrix_set (mW , i, j, 0);
			}
		}
	}

	/* Lambda  inverse */
	gsl_matrix * mLambdaInv = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots);

	int s;
	gsl_linalg_LU_decomp (mLambda, mPermutation, &s);
	gsl_linalg_LU_invert (mLambda, mPermutation, mLambdaInv);


	/* W matrix */
	gsl_matrix * mLambdaInvTheta = gsl_matrix_alloc (iNumOfKnots, iNumOfKnots);
	for(int i=0; i<iNumOfKnots; i++)
	{
		for(int j=0; j<iNumOfKnots; j++)
		{
			gsl_matrix_set ( mLambdaInvTheta  ,1,j,0 );
		}
	}

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mLambdaInv,mTheta,0.0, mLambdaInvTheta);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mP,mLambdaInvTheta,1.0, mW);

	/* WStar matrix */
	cout <<"calculate star started "<<endl;
	CalculateWStarOut( iNumOfObs, vTimes, iNumOfKnots,  vKnots, mLambdaInvTheta, mWStar );
	cout <<"calculate star finished "<<endl;

	/* WTilde matrix */
	for(int i=0; i<iNumOfObs; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			gsl_matrix_set (mWTilde ,i,j, gsl_matrix_get (mW ,i,j) - gsl_matrix_get (mW ,i,iNumOfKnots-1)*gsl_matrix_get (mWStar ,0,j)/gsl_matrix_get (mWStar ,0,iNumOfKnots-1) );
		}
	}

	gsl_matrix_free (mTheta);
	gsl_matrix_free (mLambda);
	gsl_matrix_free (mP);
	gsl_matrix_free (mW);
	gsl_matrix_free (mWStar);
	gsl_matrix_free (mLambdaInv);

	gsl_permutation_free (mPermutation);

	delete dLambda;
	delete dH;
	delete dHPlus;

}



void CalculateInitialLogVol(double* vData, double * vTime,  int iNumOfObs,  int iWindow,double * vInitialLogVol)
{
	for(int i=iWindow-1; i<iNumOfObs;i++)
	{
		double dVar=CalcVol(vData, i, iWindow);
		if(dVar==0)
		{
			vInitialLogVol[i]=log(0.01);
		}
		else
		{
			vInitialLogVol[i]=log(dVar/2);
		}
	}

	for(int i=0; i<iWindow-1;i++)
	{
		vInitialLogVol[i]=vInitialLogVol[iWindow-1]-(iWindow-1-i)*(vInitialLogVol[2*iWindow-1]-vInitialLogVol[iWindow-1])/iWindow;
	}
}


void CalculateInitialState(struct AllParam sParam, double * vInitialLogVol, int iNumOfObs, int iNumOfKnots,   gsl_matrix * mWTilde  )
{
	gsl_matrix * mYBar=gsl_matrix_alloc (iNumOfObs, 1);
	gsl_matrix * mWBar=gsl_matrix_alloc(iNumOfObs,iNumOfKnots-1);

	double * dY=new double;

	double dMeanLogInt=0;
	for(int i=0; i<iNumOfObs; i++)
	{
		dMeanLogInt=dMeanLogInt+vInitialLogVol[i]/iNumOfObs;
	}

	//int iTempCount=0;
	for(int i=0; i<iNumOfObs; i++)
	{
//			cout << "dMeanLogInt "<<dMeanLogInt << endl;
			dY[0]=vInitialLogVol[i]-dMeanLogInt;

			gsl_matrix_set (mYBar,i,0, dY[0] );
//			cout <<"mYBar "<<gsl_matrix_get (mYBar,i,0)<<endl;
			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,i,j, gsl_matrix_get (mWTilde, i, j));
//				cout <<"mWBar "<<gsl_matrix_get (mWBar,i,j)<<endl;
			}
	}

	/* Variance */
	gsl_matrix * mWW=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVar=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_matrix * mVarInv=gsl_matrix_alloc (iNumOfKnots-1, iNumOfKnots-1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mWBar,0.0, mWW);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			if(j==i)
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j)+sParam.dPriorSeasonalVar[0]);
			}
			else
			{
				gsl_matrix_set (mVarInv,i,j,gsl_matrix_get(mWW,i,j));
			}
		}
	}

	gsl_permutation * mPermutation = gsl_permutation_alloc (iNumOfKnots-1);
	int s;
	gsl_linalg_LU_decomp (mVarInv, mPermutation, &s);
	gsl_linalg_LU_invert (mVarInv, mPermutation, mVar);


	/* Mean */
	gsl_matrix * mWY=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMean=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_matrix * mMeanTemp=gsl_matrix_alloc (iNumOfKnots-1, 1);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mWBar,mYBar,0.0, mWY);

	for(int i=0; i<iNumOfKnots-1; i++)
	{
		gsl_matrix_set (mMeanTemp,i,0,gsl_matrix_get(mWY,i,0)+sParam.dPriorSeasonalMean[0]/sParam.dPriorSeasonalVar[0]);
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mVar,mMeanTemp,0.0, mMean);


	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(mMean,j,0)*gsl_matrix_get(mWTilde,i,j);
		}

	}

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vX[i]=vInitialLogVol[i]-dMeanLogInt-sParam.vS[i];
	}

	sParam.dMu[0]=dMeanLogInt;

	gsl_matrix_free (mYBar);
	gsl_matrix_free (mWBar);
	gsl_matrix_free (mWW);
	gsl_matrix_free (mWY);
	gsl_matrix_free (mMean);
	gsl_matrix_free (mMeanTemp);
	gsl_matrix_free (mVar);
	gsl_matrix_free (mVarInv);

	delete dY;

}
