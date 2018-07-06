% vTempRow.push_back(sParam.dMu[0]);
% vTempRow.push_back(sParam.dPhi[0]);
% vTempRow.push_back(sParam.dSigma2[0]);
% vTempRow.push_back(sParam.dGamma[0]);
% for(int k=0; k<iNumOfKnots-1; k++)
% {
%     vTempRow.push_back(sParam.vBeta[k]);
% }
% vTempRow.push_back(sParam.dNu[0]);


Chain = csvread('Data_29_OrdT_EstimationResults.csv');
 
param = {'mu','phi','sigma2','gamma','beta1','beta2','nu'};
theta_true = [-1.7 0.97 0.02 0.001 NaN NaN 15];

C=sum(diff(Chain5)~=0)/size(Chain5,1); % acceptance rates
%     0.9999    0.0959    0.0959    0.3841    1.0000    1.0000    0.6200


nu = sort(Chain3(:,end));
nu_u = unique(nu);

dPriorA = 9;
dPriorB = 1.5;

prior_nu = @(nu) exp(dPriorA*log(dPriorB) + dPriorA*log(nu) - dPriorB*nu - gammaln(dPriorA));
prior_nu2 = @(nu) (dPriorB^dPriorA).*(nu.^dPriorA)./((exp(nu).^dPriorB).*gamma(dPriorA));
dNu = 0.1:0.02:128;
figure(13)
plot(dNu,prior_nu(dNu));
hold on
plot(dNu,prior_nu2(dNu),'r');
hold off



figure(1);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(2,4,ii)
    hold on
    
%     plot(Chain1(:,ii),'b')
    plot(Chain(:,ii),'m')
 %     plot(Chain_old(:,ii),'k')

%     plot(theta_true(ii)+0*Chain2(:,ii),'r')
    hold off
%     title(param{ii})
end

figure(21)
prior_gamma = @(xx) betapdf(xx,1.7,10);
plot(0.0:0.01:0.5,prior_gamma(0.0:0.01:0.5))
hold on
line([0.001,0.001],[0,5],'Color','r','LineWidth',2)

figure(4)
subplot(2,2,1)
scatter(Chain2(:,2),Chain2(:,3))
subplot(2,2,2)
plot(Chain2(:,2))
title(param{2})
subplot(2,2,4)
plot(Chain2(:,3))
title(param{3})


figure(11);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(3,7,ii)
    autocorr(Chain1(:,ii),100); xlabel(''); ylabel('') ;
    title([param{ii}, 'no int div'])

    subplot(3,7,7+ii)
    autocorr(Chain5(:,ii),100); xlabel(''); ylabel('') ;
    title([param{ii}, 'split'])


    subplot(3,7,14+ii)
    autocorr(Chain_old(:,ii),100); xlabel(''); ylabel('') ;
    title([param{ii}, 'old'])
end


figure(12);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(3,7,ii)
    plot(Chain1(1:40000,ii)); xlabel(''); ylabel('') ;
    title([param{ii}, 'no int div'])

    subplot(3,7,7+ii)
    plot(Chain5(:,ii)); xlabel(''); ylabel('') ;
    title([param{ii}, 'split'])


    subplot(3,7,14+ii)
    plot(Chain_old(1:40000,ii)); xlabel(''); ylabel('') ;
    title([param{ii}, 'old'])
end


set(gcf,'PaperPositionMode','auto');
print('../src/DNB_draws_trace2.png','-dpng','-r0')


LogitTransform = @(x) log(x./(1-x));
Gamma_logit = LogitTransform(Chain1(:,4));
plot(Gamma_logit )
sum(Gamma_logit(1:60000))


ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(4,4,ii)
  
    autocorr(Chain1(:,ii),100)
    xlabel('')
    ylabel('')
    title(param{ii})
    subplot(4,4,8+ii)

    autocorr(Chain2(:,ii),100)
    xlabel('')
    ylabel('')    
    title(param{ii})

end
set(gcf,'PaperPositionMode','auto');
print(ff,'../src/DNB_draws_acf.png','-dpng','-r0')



mu = Chain(:,1);
phi = Chain(:,2);
sigma2 = Chain(:,3);
gamma = Chain(:,4);
beta1 = Chain(:,5);
beta2 = Chain(:,6);
nu = Chain(:,7);


DNB = csvread('../src/OutPut/__DNB_SimulatedData.csv');
y = DNB(:,1);
DNB_Time = DNB(:,2);


DNB2 = csvread('src/OutPut/__DNB_SimulatedLogInt.csv');
DNB_Season = csvread('src/OutPut/__DNB_vSeason.csv');
DNB_Time = csvread('src/OutPut/__DNB_vTime.csv');


S = reshape(DNB(:,2),2000,10);
S2 = reshape(DNB_Season,2000,10);


sum(DNB_Time==DNB(:,2))


DNB_short = csvread('src/DNB_Output_short/SimulatedData.csv');
6.5*3600









XHigh = csvread('../src/V3/OutPut/__DNB_vXHigh.csv');
hold on
plot(XHigh(:,ii),'k')
for ii = 1:5
    plot(XHigh(:,ii))
end
hold off