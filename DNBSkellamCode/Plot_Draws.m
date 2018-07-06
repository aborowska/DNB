SK_DRAWS = csvread('C:\Users\aga\Dropbox\New Projects\Integer\Codes\DNBSkellamCode\src\__Sk_EstimationResults.csv');
DNB_DRAWS = csvread('C:\Users\aga\Dropbox\New Projects\Integer\Codes\DNBSkellamCode\src\__DNB_EstimationResults.csv');

% Sk_Param_init = [0,0.95,0.01,0.3, NaN, NaN];
% DNB_Param_init = [0,0.95,0.01,0.3,NaN, NaN,20];
Sk_Param_init = [0,0.95,0.01,0.3, 0, 0];
DNB_Param_init = [0,0.95,0.01,0.3,0, 0,20];

Sk_Param_true = [-1.7, 0.97, 0.02, 0.001, NaN, NaN];
DNB_Param_true = [-1.7, 0.97, 0.02, 0.001, NaN, NaN, 15];

figure(1)
set(gcf, 'Units', 'normalized', 'Position', [0.1 0 0.9 1]);
for ii = 1:6
    subplot(2,3,ii)
    plot(SK_DRAWS(:,ii))
    hold on
    plot(Sk_Param_init(ii) + 0*SK_DRAWS(:,ii),'g')
    plot(Sk_Param_true(ii) + 0*SK_DRAWS(:,ii),'r')

    hold off
    title(num2str(ii))
end
% subtitle('SK')

figure(2)
set(gcf, 'Units', 'normalized', 'Position', [0.1 0 0.9 1]);
for ii = 1:7
    subplot(2,4,ii)
    plot(DNB_DRAWS(:,ii))
    hold on
    plot(DNB_Param_init(ii) + 0*DNB_DRAWS(:,ii),'g')
    plot(DNB_Param_true(ii) + 0*DNB_DRAWS(:,ii),'r')

    hold off
    title(num2str(ii))
end
% subtitle('SK')

% initial 
% 	sParam.dMu[0]=
%[0;0.95;0.01;0.3]