void DrawSesonal(struct AllParam sParam,  int iNumOfObs, int iNumOfKnots,   gsl_matrix * mWTilde ,const gsl_rng * gsl_random_num)
{

	/* Number of observations in the OLS */
	int iNumOfOLSObs=0;
	for(int i=0; i<iNumOfObs;i++)
	{
		if(sParam.vN[i]==0)
		{
			iNumOfOLSObs=iNumOfOLSObs+1;
		}
		else
		{
			iNumOfOLSObs=iNumOfOLSObs+2;
		}
	}
	gsl_matrix * mYBar=gsl_matrix_alloc (iNumOfOLSObs, 1);
	gsl_matrix * mWBar=gsl_matrix_alloc(iNumOfOLSObs,iNumOfKnots-1);

	double * dY1=new double;
	double * dY2=new double;

	double * dS1=new double;
	double * dS2=new double;


	int iTempCount=0;
	for(int i=0; i<iNumOfObs; i++)
	{


		if(sParam.vN[i]==0)
		{

			dY1[0]=sParam.mAuxY[i*2]-(sParam.dMu[0]+sParam.vX[i]);
			dS1[0]=pow(sParam.mAuxH[i*4],0.5);

			gsl_matrix_set (mYBar,iTempCount,0, dY1[0]/dS1[0] );

			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,iTempCount,j, gsl_matrix_get (mWTilde, i, j)/dS1[0]);
			}

			iTempCount=iTempCount+1;


		}
		else
		{

			dY1[0]=sParam.mAuxY[i*2]-(sParam.dMu[0]+sParam.vX[i]);
			dY2[0]=sParam.mAuxY[i*2+1]-(sParam.dMu[0]+sParam.vX[i]);
			dS1[0]=pow(sParam.mAuxH[i*4],0.5);
			dS2[0]=pow(sParam.mAuxH[i*4+3],0.5);

			gsl_matrix_set (mYBar,iTempCount,0, dY1[0]/dS1[0]);
			gsl_matrix_set (mYBar,iTempCount+1,0, dY2[0]/dS2[0]);

			for(int j=0; j<iNumOfKnots-1; j++)
			{
				gsl_matrix_set (mWBar,iTempCount,j, gsl_matrix_get (mWTilde, i, j)/dS1[0]);
				gsl_matrix_set (mWBar,iTempCount+1,j, gsl_matrix_get (mWTilde, i, j)/dS2[0]);
			}

			iTempCount=iTempCount+2;
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


	/* Draw Normal */
	gsl_linalg_cholesky_decomp (mVar);
	gsl_matrix * vNormal=gsl_matrix_alloc (iNumOfKnots-1, 1);


	string sFile;
	for(int i=0; i<iNumOfKnots-1;i++)
	{

		if(!isthisfinite(gsl_matrix_get(mMean,i,0)))
		{
			double * YBar=new double[iNumOfOLSObs];
			for(int k=0; k<iNumOfOLSObs;k++)
			{
				YBar[k]=gsl_matrix_get(mYBar,k,0);
			}
			sFile="vYBarSpineError.csv";
			WriteOutDoubleArray( YBar , iNumOfOLSObs, 1, sFile);

			double * WBar=new double[iNumOfOLSObs*(iNumOfKnots-1)];
			for(int k=0; k<iNumOfOLSObs;k++)
			{
				for(int l=0; l<iNumOfKnots-1;l++)
				{
					WBar[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mWBar,k,l);
				}
			}
			sFile="vWBarSpineError.csv";
			WriteOutDoubleArray( WBar , iNumOfOLSObs, iNumOfKnots-1, sFile);

			double * Mean=new double[iNumOfKnots-1];
			for(int k=0; k<iNumOfKnots-1;k++)
			{
				Mean[k]=gsl_matrix_get(mMean,k,0);
			}
			sFile="vMeanSpineError.csv";
			WriteOutDoubleArray( Mean , iNumOfKnots-1, 1, sFile);





			sFile="vHSpineError.csv";
			WriteOutDoubleArray( sParam.mAuxH , iNumOfObs, 4, sFile);


			sFile="vAuxYSpineError.csv";
			WriteOutDoubleArray( sParam.mAuxY , iNumOfObs, 2, sFile);

			sFile="vXSpineError.csv";
			WriteOutDoubleArray( sParam.vX , iNumOfObs, 1, sFile);

			sFile="vSSpineError.csv";
			WriteOutDoubleArray( sParam.vS , iNumOfObs, 1, sFile);

			double * Var=new double[(iNumOfKnots-1)*(iNumOfKnots-1)];
			for(int k=0; k<iNumOfKnots-1;k++)
			{
				for(int l=0; l<iNumOfKnots-1;l++)
				{
					Var[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mVar,k,l);
				}
			}
			sFile="vVarSpineError.csv";
			WriteOutDoubleArray( Var , iNumOfKnots-1, iNumOfKnots-1, sFile);

			delete [] YBar;
			delete [] WBar;
			delete [] Mean;
			delete [] Var;

			exit(1);
		}

		for(int j=0; j<iNumOfKnots-1;j++)
		{
			if(!isthisfinite(gsl_matrix_get(mVar,i,j)))
			{
				double * YBar=new double[iNumOfOLSObs];
				for(int k=0; k<iNumOfOLSObs;k++)
				{
					YBar[k]=gsl_matrix_get(mYBar,k,0);
				}
				sFile="vYBarSpineError.csv";
				WriteOutDoubleArray( YBar , iNumOfOLSObs, 1, sFile);

				double * WBar=new double[iNumOfOLSObs*(iNumOfKnots-1)];
				for(int k=0; k<iNumOfOLSObs;k++)
				{
					for(int l=0; l<iNumOfKnots-1;l++)
					{
						WBar[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mWBar,k,l);
					}
				}
				sFile="vWBarSpineError.csv";
				WriteOutDoubleArray( WBar , iNumOfOLSObs, iNumOfKnots-1, sFile);

				double * Mean=new double[iNumOfKnots-1];
				for(int k=0; k<iNumOfKnots-1;k++)
				{
					Mean[k]=gsl_matrix_get(mMean,k,0);
				}
				sFile="vMeanSpineError.csv";
				WriteOutDoubleArray( Mean , iNumOfKnots-1, 1, sFile);

				sFile="vHSpineError.csv";
				WriteOutDoubleArray( sParam.mAuxH , iNumOfObs, 4, sFile);

				sFile="vAuxYSpineError.csv";
				WriteOutDoubleArray( sParam.mAuxY , iNumOfObs, 2, sFile);

				sFile="vXSpineError.csv";
				WriteOutDoubleArray( sParam.vX , iNumOfObs, 1, sFile);

				sFile="vSSpineError.csv";
				WriteOutDoubleArray( sParam.vS , iNumOfObs, 1, sFile);

				double * Var=new double[(iNumOfKnots-1)*(iNumOfKnots-1)];
				for(int k=0; k<iNumOfKnots-1;k++)
				{
					for(int l=0; l<iNumOfKnots-1;l++)
					{
						Var[k*(iNumOfKnots-1)+l]=gsl_matrix_get(mVar,k,l);
					}
				}
				sFile="vVarSpineError.csv";
				WriteOutDoubleArray( Var , iNumOfKnots-1, iNumOfKnots-1, sFile);

				delete [] YBar;
				delete [] WBar;
				delete [] Mean;
				delete [] Var;
				exit(1);
			}
		}

	}

	/* Draw Beta */
	gsl_matrix * vBeta=gsl_matrix_alloc (iNumOfKnots-1, 1);
	for(int i=0; i<iNumOfKnots-1;i++)
	{
		gsl_matrix_set (vNormal,i,0,gsl_ran_gaussian(gsl_random_num,1));

		gsl_matrix_set (vBeta,i,0,0);
		for(int j=0; j<=i; j++)
		{
			gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(vNormal,j,0)*gsl_matrix_get(mVar,i,j));
		}

		gsl_matrix_set (vBeta,i,0, gsl_matrix_get(vBeta,i,0)+gsl_matrix_get(mMean,i,0));

		sParam.vBeta[i]=gsl_matrix_get(vBeta,i,0);
	}

	/* Seasonal */

	for(int i=0; i<iNumOfObs; i++)
	{
		sParam.vS[i]=0;
		for(int j=0; j<iNumOfKnots-1; j++)
		{
			sParam.vS[i]=sParam.vS[i]+gsl_matrix_get(vBeta,j,0)*gsl_matrix_get(mWTilde,i,j);
		}

	}



	/*print*/
//	for(int j=0; j<iNumOfKnots-1; j++)
//	{
//		cout<<"mean temp at "<<j<<" is "<<gsl_matrix_get(mMeanTemp,j,0)<<endl;
//	}
//	for(int j=0; j<iNumOfKnots-1; j++)
//	{
//		cout<<"mean at "<<j<<" is "<<gsl_matrix_get(mMean,j,0)<<endl;
//	}
//
//
//	for(int j=0; j<iNumOfKnots-1; j++)
//	{
//		cout<<"beta at "<<j<<" is "<<gsl_matrix_get(vBeta,j,0)<<endl;
//	}
//
//
//
//	for(int i=4764; i<4768;  i++)
//	{
//
//		printf ("vSesonal(%d) = %g\n", i, sParam.vS[i]);
//
//	}




	gsl_matrix_free (mYBar);
	gsl_matrix_free (mWBar);
	gsl_matrix_free (mWW);
	gsl_matrix_free (mWY);
	gsl_matrix_free (mMean);
	gsl_matrix_free (mMeanTemp);
	gsl_matrix_free (mVar);
	gsl_matrix_free (mVarInv);

	delete dS1;
	delete dS2;
	delete dY2;
	delete dY1;
}

