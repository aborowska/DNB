function CaculateDNB(vData, iNumOfObs, sParam, vDNB, vZeroDNBDens, vIndicator)
	dLambda= 0;
	dLambdaRatio=0;
	dNuRatio=0;
	iAbsX=0;

	for ii = 1:iNumOfObs	
%		cout <<"log lambda "<< <<endl;
		dLambda = exp(sParam.dMu(1)+sParam.vS(ii)+sParam.vX(ii));
		dLambdaRatio = dLambda/(dLambda + sParam.dNu(1));
		%dNuRatio(1)=sParam.dNu(1)/(dLambda(1)+sParam.dNu(1));
		dNuRatio = 1-dLambdaRatio;
		iAbsX = abs(vData(ii));

% %		cout << "Hyper original Star "<<endl;
% %		cout << "Hyper "<<gsl_sf_hyperg_2F1(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , pow(dLambdaRatio(1),2))<<endl;
% %		cout << "Hyper original End "<<endl;
% %		cout << "Hyper Star "<<endl;
% %		cout << "Hyper "<<gsl_sf_hyperg_2F1_renorm(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , pow(dLambdaRatio(1),2))<<endl;
% %		cout << "Hyper End "<<endl;
% 
% %		cout <<"log DNB "<<2*sParam.dNu(1)*log(dNuRatio(1))+iAbsX(1)*log(dLambdaRatio(1))+lgamma(sParam.dNu(1)+iAbsX(1))-lgamma(sParam.dNu(1))
% %						+log(gsl_sf_hyperg_2F1_renorm(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , pow(dLambdaRatio(1),2)))<<endl;
% 
% %		vDNB(ii)=exp(2*sParam.dNu(1)*log(dNuRatio(1))+iAbsX(1)*log(dLambdaRatio(1))+lgamma(sParam.dNu(1)+iAbsX(1))-lgamma(sParam.dNu(1))
% %				+log(gsl_sf_hyperg_2F1_renorm(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , pow(dLambdaRatio(1),2))));
% 
% %		cout <<"vData(ii) "<< vData(ii) <<endl;
% %		cout <<"Log Lambda "<<sParam.dMu(1)+sParam.vX(ii)<<endl;
% 		gsl_sf_result result;
        
        result = gsl_sf_hyperg_2F1_e(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , dLambdaRatio^2);
%These routines compute the Gauss hypergeometric function 2F1(a,b,c,x) = F(a,b,c,x) for |x| < 1. 
 

        vDNB(ii) = exp(2*sParam.dNu(1)*log(dNuRatio(1)) + ...
            iAbsX(1)*log(dLambdaRatio(1)) + ...
            lgamma(sParam.dNu(1)+iAbsX(1)) - ...
            lgamma(sParam.dNu(1)) - ...
            lgamma(iAbsX(1)+1)...
			+ log(result.val));
       

%		vDNB(ii)=exp(2*sParam.dNu(1)*log(dNuRatio(1))+iAbsX(1)*log(dLambdaRatio(1))+lgamma(sParam.dNu(1)+iAbsX(1))-lgamma(sParam.dNu(1))-lgamma(iAbsX(1)+1)
%					+log(gsl_sf_hyperg_2F1(sParam.dNu(1)+iAbsX(1),sParam.dNu(1), iAbsX(1)+1 , pow(dLambdaRatio(1),2))));

		if (vData(ii) == 0)
			vIndicator(ii) = 1;		
        else
			vIndicator(ii) = 0;
        end

%		cout << "vIndicator "<<vIndicator(ii)<<endl;
		vZeroDNBDens(ii) = sParam.dGamma(1)*vIndicator(ii) + (1-sParam.dGamma(1))*vDNB(ii) ;
	
    end
end