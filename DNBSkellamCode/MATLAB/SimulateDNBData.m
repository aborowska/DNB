% /*
%  * Simulate  from dynamic zero inflated negative binomial difference model
%  *
%  * iNumOfSim	number of simulated data points
%  *
%  * Note that the algorithm uses :
%  * dMuTrue;
%  * dPhiTrue;
%  * Sigma2True;
%  * dGammaTrue;
%  * dNuTrue;
%  *
%  */

function [vData,vTimes,vXTrue,vSeasonTrue] = ...
    SimulateDNBData(theta_true, iNumOfSim)

% 	string sType;
% 	sType="DNB";
% 	const gsl_rng_type * T;
% 	gsl_rng * r;

    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s); 
    
% 	T=gsl_rng_default;
% 	r=gsl_rng_alloc(T);
% 	gsl_rng_set(r, 1);

	vData = zeros(iNumOfSim,1);
	vXTrue  = zeros(iNumOfSim,1);

% 	iNumOfObs = iNumOfSim;

    dMuTrue = theta_true(1);
    dPhiTrue = theta_true(2);
    dSigma2True = theta_true(3);
    dGammaTrue = theta_true(4);
    dNuTrue = theta_true(5);
% 	cout<<"dMuTrue "<< dMuTrue<<endl;
% 	cout<<"dPhiTrue "<< dPhiTrue<<endl;
% 	cout<<"dSigma2True "<< dSigma2True <<endl;
% 	cout<<"dGammaTrue "<< dGammaTrue <<endl;
% 	cout<<"dNuTrue "<< dNuTrue<<endl;
	dStateSigma = sqrt(dSigma2True);

	iNumOfSimPerDay = 2000;
% 	// SimulateSeasonal(iNumOfSimPerDay, iNumOfSim/iNumOfSimPerDay,  vTimes,  vSeasonTrue);
    [vTimes, vSeasonTrue] = SimulateSeasonal_str(iNumOfSimPerDay, iNumOfSim);

    for ii = 1:iNumOfSim
        if (ii==1)		
        % 			/* Log intensity */
        % 			vXTrue(ii)= gsl_ran_gaussian(r,sqrt(dSigma2True/(1-dPhiTrue*dPhiTrue)) );		
            vXTrue(ii)= sqrt(dSigma2True/(1-dPhiTrue*dPhiTrue))*randn;		
        else		
        % 			vXTrue(ii) = dPhiTrue * vXTrue(ii-1)+ gsl_ran_gaussian(r,dStateSigma);
            vXTrue(ii) = dPhiTrue * vXTrue(ii-1) + dStateSigma * randn;
        end
% 		/* Conditional DNB return */
% 		dPoissionMean1= exp(dMuTrue +vSeasonTrue(ii) + vXTrue(ii)) * gsl_ran_gamma(r, dNuTrue, 1/dNuTrue);
% 		dPoissionMean2= exp(dMuTrue +vSeasonTrue(ii) + vXTrue(ii)) * gsl_ran_gamma(r, dNuTrue, 1/dNuTrue);
		dPoissionMean1= exp(dMuTrue +vSeasonTrue(ii) + vXTrue(ii)) * gamrnd(dNuTrue, 1/dNuTrue);
		dPoissionMean2= exp(dMuTrue +vSeasonTrue(ii) + vXTrue(ii)) * gamrnd(dNuTrue, 1/dNuTrue);
		% int IDNB
%         iDNB = gsl_ran_poisson(r,dPoissionMean1) -  gsl_ran_poisson(r,dPoissionMean2) ;
        iDNB = poissrnd(dPoissionMean1) -  poissrnd(dPoissionMean2) ;
% 		dU=gsl_rng_uniform(r);
		dU = rand;
        if (dU <= dGammaTrue)	
            vData(ii)=0;
        else	
            vData(ii) = iDNB;	
        end
		fprintf('Data at %i is %6.4f\n', ii, vData(ii));
    end

end
