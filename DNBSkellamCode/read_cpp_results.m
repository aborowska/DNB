param = {'mu','phi','sigma2','gamma','beta1','beta2','nu'};

theta_true = [-1.7 0.97 0.02 0.1 NaN NaN 15];

Chain_sim1 = csvread('src/Simulations/__DNB1_EstimationResults.csv');
Chain_sim2 = csvread('src/Simulations/__DNB2_EstimationResults.csv');
Chain_sim3 = csvread('src/Simulations/__DNB3_EstimationResults.csv');


Chain_sim4 = csvread('src/Simulations/__DNB4_EstimationResults.csv');
Chain_sim5 = csvread('src/Simulations/__DNB5_EstimationResults.csv');
Chain_sim6 = csvread('src/Simulations/__DNB6_EstimationResults.csv');
Chain_sim7 = csvread('src/Simulations/__DNB7_EstimationResults.csv');

Chain_sim6 = csvread('src/Simulations/__DNB6_EstimationResults.csv');
Chain_sim7 = csvread('src/Simulations/__DNB7_EstimationResults.csv');
Chain_sim8 = csvread('src/Simulations/__DNB8_EstimationResults.csv');
Chain_sim9 = csvread('src/Simulations/__DNB9_EstimationResults.csv');


Chain_sim10 = csvread('src/Simulations/__DNB10_EstimationResults.csv');
Chain_sim11 = csvread('src/Simulations/__DNB11_EstimationResults.csv');
Chain_sim12 = csvread('src/Simulations/__DNB12_EstimationResults.csv');
Chain_sim13 = csvread('src/Simulations/__DNB13_EstimationResults.csv');


Chain_sim14 = csvread('src/Simulations/__DNB14_EstimationResults.csv');
Chain_sim15 = csvread('src/Simulations/__DNB15_EstimationResults.csv');
Chain_sim16 = csvread('src/Simulations/__DNB16_EstimationResults.csv');
Chain_sim17 = csvread('src/Simulations/__DNB17_EstimationResults.csv');



Chain_Sk1 = csvread('src/Simulations/__Sk1_EstimationResults.csv');
figure(111);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:6
    subplot(2,4,ii)
    hold on
    plot(Chain_Sk1(1:end,ii))
    plot(theta_true(ii)+0*Chain_GT1(1:end,ii),'r')
    hold off
    title(param{ii})
end










Chain_notrans = csvread('src/V6/OutPut_notrans1/__DNB_EstimationResults.csv');
Chain_notrans_mix = csvread('src/V6/OutPut/__DNB_EstimationResults.csv');






theta_true = [-1.7 0.97 0.02 0.1 NaN NaN 10];

Chain_GT1= csvread('src/NewGT/__DNB1_EstimationResults.csv');
Chain_GTprior= csvread('src/NewGT/__DNB1_EstimationResults_prior1.csv');
Chain_GT2= csvread('src/NewGT/__DNB2_EstimationResults.csv');
Chain_GT3= csvread('src/NewGT/__DNB3_EstimationResults.csv');
 


figure(112);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(2,4,ii)
    hold on
    plot(Chain_GT1(1:end,ii),'b') 
    plot(Chain_GTprior(1:end,ii),'c') 
    plot(Chain_GT2(1:end,ii),'g') 
    plot(Chain_GT3(1:end,ii),'m') 
    plot(theta_true(ii)+0*Chain_GT1(1:end,ii),'r')
    hold off
    title(param{ii})
end


%% DEBUG skellam

Chain_Sk_001 = csvread('src/Simulations/Debug/001___EstimationResults.csv');
Chain_Sk_010 = csvread('src/Simulations/Debug/010___EstimationResults.csv');
Chain_Sk_050 = csvread('src/Simulations/Debug/050___EstimationResults.csv');

for ii=1:6
    subplot(2,4,ii)
    hold all
    plot(Chain_Sk_001(1:end,ii))
    plot(Chain_Sk_010(1:end,ii))
    plot(Chain_Sk_050(1:end,ii))
    plot(theta_true(ii)+0*Chain_Sk_001(1:end,ii),'r')
    hold off
    title(param{ii})
end






%%


Chain2 = csvread('../src/V3/OutPut/__DNB_EstimationResults.csv');

Chain5 = csvread('../src/V5/OutPut/__DNB_EstimationResults.csv');

Chain1 = csvread('src/OutPut/__DNB_EstimationResults.csv');
Chain_old = csvread('src/OutPut_2.1/__DNB_EstimationResults.csv');
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



figure(2);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(2,4,ii)
    hold on
    
%     plot(Chain1(:,ii),'b')
%     plot(Chain_notrans(:,ii),'m')
%     plot(Chain_notrans_mix(:,ii),'b')
%     plot(Chain_sim1(1:end,ii),'b')
%     plot(Chain_sim2(1:end,ii),'g')
%     plot(Chain_sim3(1:end,ii),'m')
    
    plot(Chain_GTprior(1:end,ii),'b')

%     plot(Chain_sim4(1:end,ii),'b')
%     plot(Chain_sim5(1:end,ii),'g')
%     plot(Chain_sim6(1:end,ii),'m')
%     plot(Chain_sim7(1:end,ii),'c')
%     plot(Chain_old(:,ii),'k')

%     plot(theta_true(ii)+0*Chain_sim7(:,ii),'r')
%     plot(theta_true(ii)+0*Chain_GT(20001:end,ii),'r')
    hold off
%     title(param{ii})
end

figure(12)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
set(gcf,'defaulttextinterpreter','latex');
for ii=1:7
    subplot(1,7,ii)
    hold on 
%     autocorr(Chain_sim6(20001:end,ii),100) 
    autocorr(Chain_GTprior(10001:end,ii),100) 
    hold off
    title(param{ii})
end

nu_unique = nu_gt(1);
for ii = 2:40000
    if nu_unique(end) ~= nu_gt(ii) 
        nu_unique = [nu_unique; nu_gt(ii)];
    end
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
%     subplot(2,4,ii)
    autocorr(Chain_notrans(:,ii),100); xlabel(''); ylabel('') ;
    title([param{ii}, 'split no trans'])


    subplot(3,7,14+ii)
    autocorr(Chain_notrans_mix(:,ii),100); xlabel(''); ylabel('') ;
    title([param{ii}, 'mix'])
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

 
DNB_short = csvread('src/DNB_Output_short/SimulatedData.csv');










XHigh = csvread('../src/V3/OutPut/__DNB_vXHigh.csv');
hold on
plot(XHigh(:,ii),'k')
for ii = 1:5
    plot(XHigh(:,ii))
end
hold off