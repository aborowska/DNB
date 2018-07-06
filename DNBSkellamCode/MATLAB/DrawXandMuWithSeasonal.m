function DrawXandMuWithSeasonal(iNumOfKnots, mWTilde, sParam, iNumOfObs)
% // 		DrawXandMuWithSeasonal(iNumOfKnots, mWTilde, sParam, iNumOfObs, gsl_random_num );


	iDimOfObs = 2;
	iDimOfStates = iNumOfKnots-1+2;
	mDraw = zeros(iDimOfStates*iNumOfObs);

% //	for(int i=0; i<iNumOfObs; i++)
% //	
% //		sParam.mAuxY[2*i]=sParam.mAuxY[2*i]-sParam.vS(ii);
% //		sParam.mAuxY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS(ii);
% //	
	SimulationSmootherWithSeasonal(iNumOfKnots,   mWTilde ,sParam, ...
        sParam.mAuxY,sParam.mAuxH, iNumOfObs,iDimOfObs,iDimOfStates,mDraw);

% //	double * mY=new double[2*iNumOfObs];
% //	double * mH=new double[4*iNumOfObs];
% //	for(int i=0; i<iNumOfObs; i++)
% //	
% //		mY[2*i]=sParam.mAuxY[2*i]-sParam.vS(ii);
% //		mY[2*i+1]=sParam.mAuxY[2*i+1]-sParam.vS(ii);
% //	
% //
% //	for(int i=0; i<iNumOfObs; i++)
% //	
% //		mH[4*i]=sParam.mAuxH[4*i];
% //		mH[4*i+1]=sParam.mAuxH[4*i+1];
% //		mH[4*i+2]=sParam.mAuxH[4*i+2];
% //		mH[4*i+3]=sParam.mAuxH[4*i+3];
% //	
% //
    for(int j=0; j<iNumOfKnots-1; j++)
        cout <<"beta " <<mDraw(jj)<<endl;
    end

    for(int i=0; i<iNumOfObs; i++)

        sParam.vS(ii)=0;
        for(int j=0; j<iNumOfKnots-1; j++)
            sParam.vS(ii)=sParam.vS(ii)+gsl_matrix_get(mWTilde,i,j)*mDraw(jj);
        end
    end

    for j=1:(iNumOfKnots-1)
        sParam.vBeta(jj)=mDraw(jj);
    end

    sParam.dMu = mDraw(iDimOfStates-2);

    for i=1:iNumOfObs
        sParam.vX(ii)=mDraw[iDimOfStates*i+iDimOfStates-1];
    end
end