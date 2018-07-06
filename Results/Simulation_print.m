clear all
close all

param = {'$\mu$','$\varphi$','$\sigma^2_{\eta}$','$\gamma$','$\beta_1$','$\beta_2$','$\nu$'};
param_true = [-1.7000, 0.9700, 0.0200, 0.1000, NaN, NaN, 10.0000 ];

path = '../DNBSkellamCode/src/NewGT/';
% path = '';
S = 31;
M = 40000;
BurnIn = 20000;
ChainsDNB = zeros(M-BurnIn,7,S);
ESS_DNB = zeros(7,S);
IF_DNB = zeros(7,S);
for ii = 19:S
    FileName = [path,'__DNB',num2str(ii),'_EstimationResults.csv'];
    Chain = csvread(FileName);
    Chain = Chain((BurnIn+1):M,:);
    ChainsDNB(:,:,ii) = Chain;
    [ESS_DNB(:,ii), IF_DNB(:,ii)] = ESS(Chain,0);
end

Means_DNB = squeeze(mean(ChainsDNB,1));
Std_DNB = squeeze(std(ChainsDNB,1));
Q025_DNB = squeeze(quantile(ChainsDNB,0.025,1));
Q975_DNB = squeeze(quantile(ChainsDNB,0.975,1));

Means_ESS = mean(ESS_DNB,2);
Std_ESS = std(ESS_DNB,1,2);

Means_IF = mean(IF_DNB,2);
Std_IF = std(IF_DNB,1,2);

param_true(5) = mean(mean(ChainsDNB(:,5,:)));
param_true(6) = mean(mean(ChainsDNB(:,6,:)));

print_table_results_ess(param_true, Means_DNB, Std_DNB, Q025_DNB, Q975_DNB, ...
    Means_IF, Std_IF);


S2 = 12;
M = 40000;
BurnIn = 0; %20000;
ChainsSk = zeros(M-BurnIn,6,S2);
ESS_Sk = zeros(6,S2);
IF_Sk = zeros(6,S2);
for ii = 1:S2
    FileName = [path,'__Sk',num2str(ii),'_EstimationResults.csv'];
    Chain = csvread(FileName);
    Chain = Chain((BurnIn+1):M,:);
    ChainsSk(:,:,ii) = Chain;
%     [ESS_Sk(:,ii), IF_Sk(:,ii)] = ESS(Chain,0);
end
print_table_results(param_true, ChainsDNB, ChainsSk)

%%
path = '../DNBSkellamCode/src/NewGT/';

S = 15;
M = 40000;
BurnIn = 20000;
ChainsDNB_GT = zeros(M-BurnIn,7,S);
ESS_DNB_GT = zeros(7,S);
IF_DNB_GT = zeros(7,S);
for ii = 1:S
    FileName = [path,'__DNB',num2str(ii),'_EstimationResults.csv'];
    Chain = csvread(FileName);
    Chain = Chain((BurnIn+1):M,:);
    ChainsDNB_GT(:,:,ii) = Chain;
    [ESS_DNB_GT(:,ii), IF_DNB_GT(:,ii)] = ESS(Chain,0);
end


param_true(5) = mean(mean(ChainsDNB_GT(:,5,:)));
param_true(6) = mean(mean(ChainsDNB_GT(:,6,:)));
print_table_results(param_true, ChainsDNB_GT, ChainsDNB_GT*0)




Means_DNB = squeeze(mean(ChainsDNB,1));
Std_DNB = squeeze(std(ChainsDNB,1));
Q025_DNB = squeeze(quantile(ChainsDNB,0.025,1));
Q975_DNB = squeeze(quantile(ChainsDNB,0.975,1));
 


Means_Sk = squeeze(mean(ChainsSk,1));
Std_Sk = squeeze(std(ChainsSk,1));
Q025_Sk = squeeze(quantile(ChainsSk,0.025,1));
Q975_Sk = squeeze(quantile(ChainsSk,0.975,1));



Means_DNB_GT = squeeze(mean(ChainsDNB_GT,1));
Std_DNB_GT = squeeze(std(ChainsDNB_GT,1));
Q025_DNB_GT = squeeze(quantile(ChainsDNB_GT,0.025,1));
Q975_DNB_GT = squeeze(quantile(ChainsDNB_GT,0.975,1));
 
mean(ESS_DNB,2)
std(ESS_DNB,1,2)
min(ESS_DNB,[],2)
max(ESS_DNB,[],2)


mean(IF_DNB,2)
std(IF_DNB,1,2)
min(IF_DNB,[],2)
max(IF_DNB,[],2)
