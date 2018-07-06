
function DrawNuInBlockDiscreteDNB(vData, sParam , iNumOfObs,  gsl_random_num, dResolution, dDelta)
% 		// DrawNuInBlockDiscreteDNB(vData, sParam ,  iNumOfObs, gsl_random_num, 0.2,3);


% 	/* Propose draw */
	dU = gsl_rng_uniform(gsl_random_num);
	iN = 2*dDelta/dResolution+1;
	dLogPosteriorNew = 0 ;
	dLogPosteriorOld = 0;
	dCum = 1.0/iN;
	dNuProposed = sParam.dNu(1)-dDelta;
	
    while (dCum<dU)
	
		dCum=dCum+1.0/iN;
		dNuProposed=dNuProposed+dResolution;
	

% 	cout<< "dNuProposed "<<dNuProposed<<endl;
% 	cout<< "dNuCurrent "<<sParam.dNu(1)<<endl;

	dNewNuLog=log(dNuProposed);
	dOldNuLog=log(sParam.dNu(1));

    if ((dNuProposed <= 128) && (dNuProposed >=0.1))
	
% 		/* Calculate acceptance probability */
		NuInBlockLogPosteriorDNB(&dNewNuLog, vData, sParam,  sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorNew);
		NuInBlockLogPosteriorDNB(&dOldNuLog, vData, sParam,  sParam.dPriorNuA, sParam.dPriorNuB,   iNumOfObs, dLogPosteriorOld);

% 		double dAccept=fmin(dLogPosteriorNew(1)- dLogPosteriorOld(1),0);

% 		/* Accept or reject*/
 		dLogU=log(gsl_rng_uniform(gsl_random_num));

% 		cout << "dAccept "<<dAccept<<" dLogU "<< dLogU<<endl;
        if(dLogU<=dAccept)
            sParam.dNu(1)=dNuProposed;
        end
        
    end
    
    end
end

function NuInBlockLogPosteriorDNB( dLogNu, vData, sParam, dPriorA, dPriorB, iNumOfObs, dLogPosterior)

    dNu=exp(dLogNu(1));
    dOldNu=sParam.dNu(1);
    sParam.dNu(1)=dNu;

    CaculateDNBLL(vData, iNumOfObs, sParam,dLogPosterior );

% 	sParam.dPriorNuA[0]=9;
% 	sParam.dPriorNuB[0]=1.5;
dPriorA = 15; %9;
dPriorB = 1.5;
dNu = 0.1 : 0.2 :128;
%     dLogPosterior= dLogPosterior + dPriorA(1)*log(dPriorB(1)) +dPriorA(1)*log(dNu) -dPriorB(1)*dNu-lgamma(dPriorA(1));
    dLogPrior =  dPriorA*log(dPriorB) + dPriorA*log(dNu) - dPriorB*dNu - gammaln(dPriorA);
    
dPrior = (dPriorB.^dPriorA) .* (dNu.^dPriorA)./((exp(dNu).^dPriorB) .* gamma(dPriorA));
dPrior2 = (dPriorB.^dPriorA) .* (dNu.^(dPriorA-1))./((exp(dNu).^dPriorB) .* gamma(dPriorA));
    sParam.dNu=dOldNu;
end

function CaculateDNBLL(const double * vData,  int iNumOfObs,  struct AllParam sParam,double * vDNBLL )

	dLambda = 0;
	dLambdaRatio = 0;
	dNuRatio = 0;
	iAbsX = 0;
	dDNB = 0;
	vDNBLL(1) = 0;

    for ii=1:iNumOfObs 
	
		dLambda(1)=exp(sParam.dMu(1)+sParam.vS(ii)+sParam.vX(ii));
		dLambdaRatio(1)=dLambda(1)/(dLambda(1)+sParam.dNu(1));
		dNuRatio(1)=1-dLambdaRatio(1);
		iAbsX(1)=abs(vData(ii));

		dDNB(1)= exp(2*sParam.dNu(1)*log(dNuRatio(1))+iAbsX(1)*log(dLambdaRatio(1))+lgamma(sParam.dNu(1)+iAbsX(1))-lgamma(sParam.dNu(1))-lgamma(iAbsX(1)+1) ...
					+log(gsl_sf_hyperg_2F1(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , dLambdaRatio(1)^2)));

        if (vData(ii)==0)
			vDNBLL(1)=vDNBLL(1)+log(sParam.dGamma(1)+(1-sParam.dGamma(1))*dDNB(1)) ;	
        else		
			vDNBLL(1)=vDNBLL(1)+log((1-sParam.dGamma(1))*dDNB(1)) ;
        end
    end
	 
end
