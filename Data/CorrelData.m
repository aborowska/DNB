start_date = '01012008';
end_date = '31122017';

tickers = {'^GSPC' 'IBM' 'JPM' 'KO'};
tickers2 = {'AAPL','F','XRX'};

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data2 = hist_stock_data(start_date, end_date, tickers2);
data_AdjClose2 = [data2([1:end]).AdjClose];

% data_Returns = 100*diff(log(data_AdjClose));
% csvwrite('Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv',data_Returns);

start_date1 = '01102008';
end_date1 = '31102008';

dataOct = hist_stock_data(start_date1, end_date1, tickers);
data_AdjCloseOct = [dataOct([1:end]).AdjClose];
 
corrcoef(data_AdjCloseOct(:,1),data_AdjCloseOct(:,2)) %0.9758
corrcoef(data_AdjCloseOct(:,1),data_AdjCloseOct(:,3)) % 0.9181
corrcoef(data_AdjCloseOct(:,1),data_AdjCloseOct(:,4)) % 0.9458

GSPC = data_AdjClose(:,1);
IBM = data_AdjClose(:,2);
JPM = data_AdjClose(:,3);
KO = data_AdjClose(:,4);

AAPL = data_AdjClose2(:,1);
F = data_AdjClose2(:,2);
XRX = data_AdjClose2(:,3);

corrcoef(GSPC,IBM) %0.5887

corrcoef(GSPC,JPM) % 0.9422

corrcoef(GSPC,KO) % 0.9460

corrcoef(GSPC,AAPL) %0.9442
 
corrcoef(GSPC,XRX) % 0.6905

corrcoef(GSPC,F) %0.6997

subplot(4,2,1)
plot(GSPC)
subplot(4,2,3)
plot(IBM)
subplot(4,2,5)
plot(JPM)
subplot(4,2,7)
plot(KO)

subplot(4,2,2)
plot(GSPC)
subplot(4,2,4)
plot(AAPL)
subplot(4,2,6)
plot(F)
subplot(4,2,8)
plot(XRX)
