% aaa = who('-regexp','..IBM2010');
% K = size(aaa,1);
Models = {'OrdN', 'OrdT', 'Sk', 'DNB'};
Data = {'2008_IBM','2008_KO','2008_JPM','2010_IBM','2010_KO','2010_JPM'};
param = {'$\mu$','$\varphi$','$\sigma^2_{\eta}$','$\gamma$','$\beta_1$','$\beta_2$','$\nu$'};

col = [      0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];

ind = zeros(4,1);
BurnIn = 20000;


jj = 1;
try
    Chain_OrdN_bid  = csvread(['../OrderedCode/src/BidPrice/',Data{jj},'_OrdNormal_EstimationResults.csv']);
    ind(1) = 1;   
catch
    Chain_OrdN_bid = NaN*ones(40000,6);
    ind(1) = 0;      
end
try
    Chain_OrdT_bid  = csvread(['../OrderedCode/src/BidPrice/',Data{jj},'_OrdT_EstimationResults.csv']);
    ind(2) = 1;
catch
    Chain_OrdT_bid = NaN*ones(40000,7);
    ind(2) = 0;      
end
try
    Chain_Sk_bid  = csvread(['../DNBSkellamCode/src/BidPrice/',Data{jj},'_Sk_EstimationResults.csv']);
    ind(3) = 1;
catch
    Chain_Sk_bid = NaN*ones(40000,6);
    ind(3) = 0;      
end
try
    Chain_DNB_bid  = csvread(['../DNBSkellamCode/src/BidPrice/',Data{jj},'_DNB_EstimationResults.csv']);
    ind(4) = 1;
catch
    Chain_DNB_bid = NaN*ones(40000,7);
    ind(4) = 0;    
end

% for ii = 1:4 
%     ind(ii) = ~isempty(who('-regexp',['..',Models{ii}]));
% end

Mean_emp = zeros(7,4);
Std_emp = zeros(7,4);
Q025_emp = zeros(7,4);
Q975_emp = zeros(7,4);
IF_emp = zeros(7,4);

for ii = 1:4
    if ind(ii)
        Chain = eval(char(who('-regexp',['..',Models{ii}])));
        if (mod(ii,2) == 1)
            Mean_emp(1:6,ii) = mean(Chain((BurnIn+1):end,:));
            Std_emp(1:6,ii) = std(Chain((BurnIn+1):end,:));
            Q025_emp(1:6,ii) = quantile(Chain((BurnIn+1):end,:),0.025);
            Q975_emp(1:6,ii) = quantile(Chain((BurnIn+1):end,:),0.975);
            [~, IF_emp(1:6,ii)] = ESS(Chain((BurnIn+1):end,:),0);
        else
            Mean_emp(1:7,ii) = mean(Chain((BurnIn+1):end,:));
            Std_emp(1:7,ii) = std(Chain((BurnIn+1):end,:));
            Q025_emp(1:7,ii) = quantile(Chain((BurnIn+1):end,:),0.025);
            Q975_emp(1:7,ii) = quantile(Chain((BurnIn+1):end,:),0.975);   
            [~, IF_emp(1:7,ii)] = ESS(Chain((BurnIn+1):end,:),0);         
        end
    end
end


print_table_results_emp(Data{jj},Mean_emp, Std_emp, Q025_emp, Q975_emp, IF_emp);

    
ff = figure(123);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.45 0.9]);
set(gcf,'defaulttextinterpreter','latex');

for ii = 1:7
    subplot(4,2,ii)
    hold on                   
    if (ii < 7)
        plot(Chain_OrdN_bid((BurnIn+1):end,ii),'color',col(2,:))          
    end
    plot(Chain_OrdT_bid((BurnIn+1):end,ii),'color',col(7,:))
    if (ii < 7)
        plot(Chain_Sk_bid((BurnIn+1):end,ii),'color',col(6,:))          
    end
    plot(Chain_DNB_bid((BurnIn+1):end,ii),'color',col(1,:))
    hold off
    title(param{ii})
    set(gca,'TickLabelInterpreter','latex');   
    if (ii == 1)
        ll = legend('OrdN','Ord$t$','Sk','$\Delta$NB');
        set(ll,'Interpreter','latex','Location','SouthEast')
    end
end


ff = figure(234);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.45 0.9]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,2,1)
autocorr(Chain_OrdN_bid((BurnIn+1):end,ii),100)   
set(gca,'TickLabelInterpreter','latex');   
title('OrdN')

subplot(2,2,2)
autocorr(Chain_OrdT_bid((BurnIn+1):end,ii),100) 
set(gca,'TickLabelInterpreter','latex');   
title('OrdT')

subplot(2,2,3)
autocorr(Chain_Sk_bid((BurnIn+1):end,ii),100) 
set(gca,'TickLabelInterpreter','latex');   
title('Sk')

subplot(2,2,4)
autocorr(Chain_DNB_bid((BurnIn+1):end,ii),100)  
set(gca,'TickLabelInterpreter','latex');   
title('DNB')


set(gcf,'PaperPositionMode','auto');
plotname = [Data{jj},'_trace.png'];
print(ff,plotname,'-dpng','-r0')
plotname = [Data{jj},'_trace.eps'];
print(ff,plotname,'-depsc','-r0')