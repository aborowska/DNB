
function DrawNuInBlockDiscreteDNB(vData, sParam , iNumOfObs,  gsl_random_num, dResolution, dDelta)
% 		// DrawNuInBlockDiscreteDNB(vData, sParam ,  iNumOfObs, gsl_random_num, 0.2, 3);


% 	/* Propose draw */
	dU = gsl_rng_uniform(gsl_random_num);
	iN = 2*dDelta/dResolution+1; % dDelta = 3; dResolution = 0.2;
	dLogPosteriorNew = 0 ;
	dLogPosteriorOld = 0;
	dCum = 1.0/iN;
	dNuProposed = sParam.dNu - dDelta;
	
    while (dCum<dU)
		dCum = dCum + 1.0/iN;
		dNuProposed = dNuProposed + dResolution;
    end
    % 	cout<< "dNuProposed "<<dNuProposed<<endl;
    % 	cout<< "dNuCurrent "<<sParam.dNu(1)<<endl;

    dNewNu =dNuProposed;
    dOldNu = sParam.dNu(1);
%         dNewNuLog=log(dNuProposed);
%         dOldNuLog=log(sParam.dNu(1));

    if ((dNuProposed <= 128) && (dNuProposed >=0.1))

% 		/* Calculate acceptance probability */
        dLogPosteriorNew = NuInBlockLogPosteriorDNB(dNewNu, vData, sParam, ...
            sParam.dPriorNuA, sParam.dPriorNuB, iNumOfObs);
        dLogPosteriorOld = NuInBlockLogPosteriorDNB(dOldNu, vData, sParam, ...
            sParam.dPriorNuA, sParam.dPriorNuB, iNumOfObs);

        dAccept=fmin(dLogPosteriorNew(1)- dLogPosteriorOld(1),0);

% 		/* Accept or reject*/
        dLogU = log(gsl_rng_uniform(gsl_random_num));

% 		cout << "dAccept "<<dAccept<<" dLogU "<< dLogU<<endl;
        if(dLogU<=dAccept)
            sParam.dNu(1)=dNuProposed;
        end

    end
end

% function dLogPosterior = NuInBlockLogPosteriorDNB(dLogNu, vData, sParam, ...
%     dPriorA, dPriorB, iNumOfObs)
function dLogPosterior = NuInBlockLogPosteriorDNB(dNu, vData, sParam, ...
    dPriorA, dPriorB, iNumOfObs)
%     dNu = exp(dLogNu);
    dOldNu = sParam.dNu;
    sParam.dNu = dNu;

    dLogPosterior = CaculateDNBLL(vData, iNumOfObs, sParam);

    dLogPosterior = dLogPosterior + dPriorA*log(dPriorB) + ...
        dPriorA*log(dNu) - dPriorB*dNu - lgamma(dPriorA);
    sParam.dNu = dOldNu;
end

function vDNBLL = CaculateDNBLL(vData, iNumOfObs, sParam)
 
% 	dLambda = 0;
% 	dLambdaRatio = 0;
% 	dNuRatio = 0;
% 	iAbsX = 0;
% 	dDNB = 0;
	vDNBLL = 0;

    for ii = 1:iNumOfObs 
	
		dLambda = exp(sParam.dMu+sParam.vS(ii)+sParam.vX(ii));
		dLambdaRatio = dLambda/(dLambda+sParam.dNu);
		dNuRatio = 1-dLambdaRatio;
		iAbsX = abs(vData(ii));
        
        % the probability mass function at 
		dDNB = exp(2*sParam.dNu*log(dNuRatio) + ...
            iAbsX*log(dLambdaRatio) + lgamma(sParam.dNu+iAbsX) - ...
            lgamma(sParam.dNu) - lgamma(iAbsX+1) ...
					+log(hypergeom([sParam.dNu+iAbsX,sParam.dNu], iAbsX+1 , dLambdaRatio^2)));
% 					+log(gsl_sf_hyperg_2F1(sParam.dNu+iAbsX,sParam.dNu, iAbsX+1 , dLambdaRatio^2)));

        if (vData(ii)==0)
			vDNBLL = vDNBLL + log(sParam.dGamma + (1-sParam.dGamma)*dDNB) ;	
        else		
			vDNBLL = vDNBLL + log((1-sParam.dGamma)*dDNB) ;
        end
    end
	 
end
