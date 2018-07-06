clear all
close all


s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

% y = load('../Data/Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv');

% [vData,vTimes,vXTrue,vSeasonTrue] = SimulateDNBData(theta_true, iNumOfSim);
% y = vData;
% h_true = vXTrue;
y = csvread('__DNB_SimulatedData.csv');
time = y(:,2);
y = y(:,1);
h_true = csvread('__DNB_SimulatedLogInt.csv');
v_DNB = csvread('__DNB_vVolEst.csv');
s_DNB = csvread('__DNB_vSeasonEst.csv');
x_DNB = csvread('__DNB_vXEst.csv');

T = length(y);

subplot(3,1,1)
plot(h_true)
subplot(3,1,2)
plot(x_DNB)
subplot(3,1,3)  
plot(v_DNB)


hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2]; 
hyper.G = [1.7, 10];

M = 10000;
BurnIn = 1000 + 3000 + 6000;
h_init = log(var(y))*ones(1,T); 

delta.t = [0.1, 0.01, 0.01, 0.01];
delta.h = 0.1;

 
h = h_init;
[theta, delta] = SV_initialise(arg0, 2);

H_DA_RW = zeros(M,T);
theta_DA_RW = zeros(M,3);

accept_DA_RW = zeros(M,1);
A_H_DA_RW = zeros(M,1);
A_theta_DA_RW = zeros(M,3);

% theta = theta_true;

tic
for ii = -BurnIn:1:M 
    if (mod(ii,1000) == 1)
        toc;
    end
    [h, acc, A_h] = update_h(y, h, theta, delta.h);
    [theta, A_theta] = update_theta_RW(h, theta, hyper, delta.t);
    if (ii > 0)
        theta_DA_RW(ii,:) = theta;
        H_DA_RW(ii,:) = h;

        accept_DA_RW(ii,1) = acc;
        A_H_DA_RW(ii,1) = A_h;
        A_theta_DA_RW(ii,:) = A_theta;        
    end
end
time_DA_RW = toc;

accept_DA_RW = accept_DA_RW/T;
A_H_DA_RW = A_H_DA_RW/T;
mean_H_DA_RW = mean(H_DA_RW);

ESS_theta_DA_RW_40 = ESS(theta_DA_RW,40);
ESS_theta_DA_RW_100 = ESS(theta_DA_RW,100);
ESS_theta_DA_RW_1000 = ESS(theta_DA_RW,1000);
ESS_theta_DA_RW_sig = ESS(theta_DA_RW,0);

ESS_H_DA_RW_40 = ESS(H_DA_RW,40);
ESS_H_DA_RW_100 = ESS(H_DA_RW,100);
ESS_H_DA_RW_1000 = ESS(H_DA_RW,1000);
ESS_H_DA_RW_sig = ESS(H_DA_RW,0);

H_subset_DA_RW = H_DA_RW(:,(2:20)*100);        

name = [path,'SV_results_param_DA_RW',data_name,'.mat'];
save(name,...
    'y','h_true','theta_true','delta','hyper',...
    'mean_H_DA_RW',...
    'theta_DA_RW',...
    'accept_DA_RW',...
    'A_H_DA_RW',...
    'A_theta_DA_RW',...
    'H_subset_DA_RW',...
    'ESS_theta_DA_RW_40','ESS_theta_DA_RW_100','ESS_theta_DA_RW_1000','ESS_theta_DA_RW_sig',...
    'ESS_H_DA_RW_40','ESS_H_DA_RW_100','ESS_H_DA_RW_1000','ESS_H_DA_RW_sig',...       
    'time_DA_RW');

