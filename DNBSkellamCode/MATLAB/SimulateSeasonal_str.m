function [vTimes,  vSeason] = ...
    SimulateSeasonal_str(iNumOfSimPerDay,  iNumOfSim)
 
    iNumOfDays = iNumOfSim/iNumOfSimPerDay;
	vSeason = zeros(iNumOfSim,1);
 	vTimes = zeros(iNumOfSim,1);
       
	iTimeInSec = 23400;
	dDuration = iTimeInSec/iNumOfSimPerDay;

	vUnnormalizedSeason = zeros(iNumOfSim,1);
% 	dDailyAvgSeason = 0;

% 	/* Unnormalized sesonal pattern */
    for d = 1:iNumOfDays
		dDailyAvgSeason=0;
        for ii = 1:iNumOfSimPerDay
			vTimes(ii+(d-1)*iNumOfSimPerDay) = ii*dDuration;
            
			vUnnormalizedSeason(ii+(d-1)*iNumOfSimPerDay) = ...
                (vTimes(ii)-10800)*(vTimes(ii)-10800)/(58320000);
            
			dDailyAvgSeason = dDailyAvgSeason + vUnnormalizedSeason(ii);
        end

		dDailyAvgSeason = dDailyAvgSeason/iNumOfSimPerDay;

% 		/* Zero mean normalized seasonal pattern */
        for ii = 1:iNumOfSimPerDay
			vSeason(ii+(d-1)*iNumOfSimPerDay) = ...
                vUnnormalizedSeason(ii+(d-1)*iNumOfSimPerDay)-dDailyAvgSeason(1);
% //			vSeason[i+d*iNumOfSimPerDay]=0;
        end
    end
end