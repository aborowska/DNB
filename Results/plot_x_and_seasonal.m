name = 'Data_2008103_9_KO_Bid';
Data = csvread([name,'.csv']);

trade = Data(:,1);
ind = zeros(1,6);
ind(1) = 1;
jj = 1;
for ii = 2:length(trade);
    if (trade(ii) < trade(ii-1))
        jj = jj+1;
        display(ii);
        ind(jj) = ii;
    end
end
ind(6) = length(trade)+1;

for ii = 1:5
    name_curr = [name,'_',num2str(ii),'.csv'];
    Data_curr = Data(ind(ii):ind(ii+1)-1,:);
    csvwrite(name_curr,Data_curr);
end

% //	vKnots[0]=0; /* 9.30 */
% //	vKnots[1]=5400; /* 11 */
% //	vKnots[2]=10800; /* 12.30 */
% //	vKnots[3]=18000;/* 14.30 */
% //	vKnots[iNumOfKnots-1]=23400; /* 16.00 */
% 
% 	vKnots[0]=300; /* 9.35 */
% //	vKnots[1]=5400; /* 11 */
% 	vKnots[1]=10800; /* 12.30 */
% //	vKnots[3]=18000;/* 14.30 */
% //	vKnots[iNumOfKnots-1]=23100; /* 15.55 */
% 	vKnots[iNumOfKnots-1]=23400; /* 16.00 */


% 10,11,12,13,14,
% 6.5*60*60 = 23400 seconds

% 5min = 5*60 = 300sec
% 30min = 30*60 = 1800 sec
% 90min = 90*60 = 2700 sec


s_KO08 = csvread('H:/Desktop/Seasonal/2008_KO_DNB_vSeasonEst.csv');
s_KO10 = csvread('H:/Desktop/Seasonal/2010_KO_DNB_vSeasonEst.csv');
s_IBM08_short = csvread('H:/Desktop/Seasonal/2008_IBM_DNB_vSeasonEst_short.csv');
s_IBM08 = csvread('H:/Desktop/Seasonal/2008_IBM_DNB_vSeasonEst.csv');
s_IBM10 = csvread('H:/Desktop/Seasonal/2010_IBM_DNB_vSeasonEst.csv');
x_KO08 = csvread('H:/Desktop/Seasonal/2008_KO_DNB_vXEst.csv');
x_KO10 = csvread('H:/Desktop/Seasonal/2010_KO_DNB_vXEst.csv');
x_IBM08_short = csvread('H:/Desktop/Seasonal/2008_IBM_DNB_vXEst_short.csv');
x_IBM08 = csvread('H:/Desktop/Seasonal/2008_IBM_DNB_vXEst.csv');
x_IBM10 = csvread('H:/Desktop/Seasonal/2010_IBM_DNB_vXEst.csv');

close all

ff = figure(1);
hold on
% subplot(2,1,1)
plot(x_IBM08)
% subplot(2,1,2)
plot(s_IBM08,'linewidth',2)
xlim([0,length(s_IBM08)])
title('IBM 2008')
hold off
print(ff,'x_s_IBM08.png','-dpng','-r0')


ff = figure(10);
hold on
% subplot(2,1,1)
plot(x_IBM08_short)
% subplot(2,1,2)
plot(s_IBM08_short,'linewidth',2)
xlim([0,length(s_IBM08_short)])
title('IBM 2008_short')
hold off
print(ff,'x_s_IBM08_short.png','-dpng','-r0')



ff = figure(2);
hold on
% subplot(2,1,1)
plot(x_IBM10)
% subplot(2,1,2)
plot(s_IBM10,'linewidth',2)
xlim([0,length(s_IBM10)])
title('IBM 2010')
hold off
print(ff,'x_s_IBM10.png','-dpng','-r0')



ff = figure(3);
hold on
% subplot(2,1,1)
plot(x_KO10)
% subplot(2,1,2)
plot(s_KO10,'linewidth',2)
xlim([0,length(s_KO10)])
title('KO 2010')
hold off
print(ff,'x_s_KO10.png','-dpng','-r0')

ff = figure(4);
hold on
% subplot(2,1,1)
plot(x_KO08)
% subplot(2,1,2)
plot(s_KO08,'linewidth',2)
xlim([0,length(s_KO08)])
title('KO 2008')
hold off
print(ff,'x_s_KO08.png','-dpng','-r0')
