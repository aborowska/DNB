function DrawGammaAdaptiveRWDNB(iNumOfObs, iNumOfIter,  gsl_random_num, sParam, vIndicator, vDNB, dCovar,  dSum)


% 	/* Proposal */
	dOmega1 = 0.05;
% 	dU=gsl_rng_uniform(gsl_random_num);
	dU = rand;
	dNewGammaLogit = 0;
	dSigma = 0;
	dLogPosteriorNew = 0;
	dLogPosteriorOld = 0;

    if (iNumOfIter<=1)
		dSigma=0.1;
    elseif (iNumOfIter==2)
		dSigma = pow(dCovar(1)-0.5*dSum(1)*dSum(1),0.5);
        if (dSigma==0)
			dSigma=0.1;
        end
    else	
		dSigma=pow(dCovar(1), 0.5);
        if (dSigma==0)
            dSigma=0.1;
        end
    end
    
    if (dU<=dOmega1)
% 		/* Static proposal */
		dNewGammaLogit = LogitTransform(sParam.dGamma(1))+gsl_ran_gaussian(gsl_random_num,0.1);	
    else	
% 		/* Proposal based on empirical covariance */
		dNewGammaLogit = LogitTransform(sParam.dGamma(1))+gsl_ran_gaussian(gsl_random_num,dSigma*2.38);
    end

% 	/* Acceptance Rate */

	 dOldGammaLogit = LogitTransform(sParam.dGamma(1));

% 	GammaLogPosteriorDNB(&dNewGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vDNB, vIndicator, dLogPosteriorNew)	;
% 	GammaLogPosteriorDNB(&dOldGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB,   iNumOfObs,vDNB, vIndicator, dLogPosteriorOld)	;
	 dLogPosteriorNew = GammaLogPosteriorDNB(dNewGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB, iNumOfObs, vDNB, vIndicator);
	 dLogPosteriorOld = GammaLogPosteriorDNB(dOldGammaLogit, sParam.dPriorGammaA, sParam.dPriorGammaB, iNumOfObs, vDNB, vIndicator);

	 dLogU = log(rand);
     dAccept = fmin(dLogPosteriorNew(1)-dLogPosteriorOld(1),0);

    if ((dLogU<=dAccept) && (isthisfinite(LogitTransformBack(dNewGammaLogit))) )
		sParam.dGamma(1)=LogitTransformBack(dNewGammaLogit);
    end

% 	/* Covariance */
    if (iNumOfIter<=1)
		dCovar(1)=dCovar(1)+pow(LogitTransform(sParam.dGamma(1)),2);
		dSum(1)=dSum(1)+LogitTransform(sParam.dGamma(1));
    elseif(iNumOfIter==2)
		dCovar(1)=(iNumOfIter-1)*pow(dSigma,2)/iNumOfIter+dSum(1)*dSum(1)/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma(1)),2)/iNumOfIter;
		dSum(1)=dSum(1)+LogitTransform(sParam.dGamma(1));
		dCovar(1)=dCovar(1)-dSum(1)*dSum(1)/(iNumOfIter*(iNumOfIter+1));	
    else	
		dCovar(1)=(iNumOfIter-1)*dCovar(1)/iNumOfIter+dSum(1)*dSum(1)/(iNumOfIter*iNumOfIter)+pow(LogitTransform(sParam.dGamma(1)),2)/iNumOfIter;
		dSum(1)=dSum(1)+LogitTransform(sParam.dGamma(1));
		dCovar(1)=dCovar(1)-dSum(1)*dSum(1)/(iNumOfIter*(iNumOfIter+1));
    end
 
end


function dLogPosterior = GammaLogPosteriorDNB(dLogGamma, dPriorA, dPriorB,  iNumOfObs, vDNB, vIndicator)

	dLogPosterior = 0;
	dGamma = LogitTransformBack(dLogGamma(1));

    for ii= 1:iNumOfObs
        dLogPosterior = dLogPosterior + log( dGamma*vIndicator(ii)+(1-dGamma)*vDNB(ii) );
    end
	dLogPosterior = dLogPosterior + (dPriorA(1) -1)*log(dGamma)+(dPriorB(1) -1)*log(1-dGamma);
end
