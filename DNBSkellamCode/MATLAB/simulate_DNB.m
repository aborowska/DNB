theta_true = [-1.7 0.97 0.02 0.001 15];
iNumOfSim = 10000;
iNumOfSimPerDay = 2000;

% vData = zeros(iNumOfSim,1);
% vTimes = zeros(iNumOfSim,1);
% vXTrue = zeros(iNumOfSim,1);
% % vNTrue = zeros(iNumOfSim,1);
% vTau1True = zeros(iNumOfSim,1);
% vTau2True = zeros(iNumOfSim,1);
% vSeasonTrue = zeros(iNumOfSim,1);

[vData,vTimes,vXTrue,vSeasonTrue] = ...
    SimulateDNBData(theta_true, iNumOfSim);


