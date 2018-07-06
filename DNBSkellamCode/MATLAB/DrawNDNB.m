function DrawNDNB(sParam, vData, iNumOfObs, vLogSkellamZ)

% 	 dCumU =  0;
% 	 dLambda =  0;
% 	 dZSum = 0;
% 	 dZ1Ratio= 0;
% 	 dZ2Ratio= 0;
% 	 dZLambda= 0;
% 	 iTempN = 0;
% 	 iNextBinom= 0;
% 	 dU= 0;


    for ii = 1:iNumOfObs
	
% //		cout<< "Draw N "<< i <<endl;
		dCumU=0;
		dLambda = exp(sParam.dMu+sParam.vS(ii)+sParam.vX(ii));
		iTempN = abs(vData(ii));
		dZSum = sParam.vZ1(ii)+sParam.vZ2(ii);
		dZ1Ratio = sParam.vZ1(ii)/dZSum;
		dZ2Ratio = sParam.vZ2(ii)/dZSum;
		dZLambda = dZSum*dLambda;
% //		cout<<"vData(ii) "<< vData(ii)<<endl;
% //		cout<<"dZSum "<<dZSum<<endl;
% //		cout<<"dZ1Ratio "<<dZ1Ratio<<endl;
% //		cout<<"dZ2Ratio "<<dZ2Ratio<<endl;
% //		cout<<"dZLambda "<<dZLambda<<endl;

        if (vData(ii)==0)
			dCumU = exp(log(sParam.dGamma)+iTempN*log(dZLambda)-dZLambda -lgamma(iTempN+1)-vLogSkellamZ(ii))+...
					exp(log(1-sParam.dGamma)+iTempN*log(dZLambda)-dZLambda+ ...
                    0.5*(iTempN+vData(ii))*log(dZ1Ratio)+0.5*(iTempN-vData(ii))*log(dZ2Ratio) ... 
					     -lgamma(0.5*(iTempN+vData(ii))+1)-lgamma(0.5*(iTempN-vData(ii))+1)-vLogSkellamZ(ii));
		
        else
			dCumU=exp(log(1-sParam.dGamma)+iTempN*log(dZLambda)- dZLambda + ...
                0.5*(iTempN+vData(ii))*log(dZ1Ratio)+0.5*(iTempN-vData(ii))*log(dZ2Ratio) ...
				     -lgamma(0.5*(iTempN+vData(ii))+1)-lgamma(0.5*(iTempN-vData(ii))+1)-vLogSkellamZ(ii));
		
        end

		iNextBinom=iTempN+2;
		dU = gsl_rng_uniform(gsl_random_num);

% //		cout << "CumU " << dCumU << " dU "<< dU<< endl;
        while (dCumU<dU)
% //		while(dCumU<1)
		
% 			/*New N */
			iTempN=iTempN+1;

            if (vData(ii)==0)		
                if(iNextBinom==iTempN )
					dCumU=dCumU+exp(log(sParam.dGamma)+iTempN*log(dZLambda)-...
                        dZLambda-lgamma(iTempN+1)-vLogSkellamZ(ii))+...
							exp(log(1-sParam.dGamma)+iTempN*log(dZLambda)-...
                            dZLambda+0.5*(iTempN+vData(ii))*log(dZ1Ratio)+...
                            0.5*(iTempN-vData(ii))*log(dZ2Ratio)...
							     -lgamma(0.5*(iTempN+vData(ii))+1)-...
                                 lgamma(0.5*(iTempN-vData(ii))+1)-vLogSkellamZ(ii));
					iNextBinom=iTempN+2;				
                else		
					dCumU=dCumU+exp(log(sParam.dGamma)+...
                        iTempN*log(dZLambda)-dZLambda-lgamma(iTempN+1)-...
                        vLogSkellamZ(ii));			
                end		
            else		
                if(iNextBinom==iTempN )				
					dCumU=dCumU+exp(log(1-sParam.dGamma)+iTempN*log(dZLambda)-...
                        dZLambda+0.5*(iTempN+vData(ii))*log(dZ1Ratio)+...
                        0.5*(iTempN-vData(ii))*log(dZ2Ratio)...
						     -lgamma(0.5*(iTempN+vData(ii))+1)-...
                             lgamma(0.5*(iTempN-vData(ii))+1)-vLogSkellamZ(ii));
					iNextBinom=iTempN+2;				
                end
            end
		

            sParam.vN(ii)=iTempN;

            if ~isfinite(sParam.vN(ii))

    % 			cout <<"N enter"<<endl;
                vNError = zeros(11,1);
                vNError(1)=ii;
                vNError(2)=sParam.vN(ii);
                vNError(3)=dCumU;
                vNError(4)=dU;
                vNError(5)=dLambda;
                vNError(6)=vLogSkellamZ(ii);
                vNError(7)=vLogSkellamZ(ii);
                vNError(8)=vData(ii);
                vNError(9)=sParam.dMu;
                vNError(10)=sParam.vS(ii);
                vNError(11)=sParam.vX(ii);
    % 			cout <<"N write"<<endl;
    % 			string sFile="vNError.csv" ;
                WriteOutDoubleArray( vNError , 11, 1, sFile);
    % 			cout <<"N del"<<endl;
            end
        end
    end
end