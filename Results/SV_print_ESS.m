clear all
EXT = {'DA','DA_RW_eff','HMM_eff_Nbin30','HMM_adapt_eff_Nbin30'};
METH = {'DA KSC', 'DA RW', 'HMM fix','HMM adapt'};
params = {'\mu','\phi','\sigma^2'};
DATA = {'', '_GSPC' , '_MSFT', '_IBM'};
data_name = {'simulated','S\\&P500','MSFT','IBM'};


mean_ESS_theta = zeros(4,4,3,4);
mean_ESS_H = zeros(4,4,4);
std_ESS_H = zeros(4,4,4);


mean_AR_theta = zeros(4,3,4);
mean_AR_H = zeros(4,4);

mean_time = zeros(4,4);

for dd = [1,2,4] % 1:4 % DATA
    data = DATA{dd};
    if ~isempty(char(data))
        path_results = 'Results/Empirical';
        path_figures = 'figures/Empirical';
    else 
        path_results = 'Results/Simulation';
        path_figures = 'figures/Simulation';
    end

    for jj = 1:4 % METHOD
        ext = EXT{jj};
        clear -regexp ^ESS_ ^time ^A_
        load([path_results,'/SV_results_param_',ext,char(data),'.mat'],...
            '-regexp','^ESS_','^time','^A_')
        aaa = who('-regexp','^ESS_theta_');
        bbb = who('-regexp','^ESS_H_');
        ccc = who('-regexp','^A_theta_');
        ddd = who('-regexp','^A_H_');        
        eee = who('-regexp','^time_'); 
        for ii = 1:4 % ESS LAG
            mean_ESS_theta(ii,jj,:,dd) = eval(char(aaa{5-ii}));
            mean_ESS_H(ii,jj,dd) = mean(eval(char(bbb{5-ii})));
            mean_ESS_H(ii,jj,dd) = mean(eval(char(bbb{5-ii})));
            std_ESS_H(ii,jj,dd) = std(eval(char(bbb{5-ii})));
        end
        if (jj > 1)
            mean_AR_theta(jj,:,dd) = mean(eval(char(ccc)));
        end
        mean_AR_H(jj,dd) = mean(eval(char(ddd)));
        mean_time(jj,dd) = mean(eval(char(eee)));       
    end  
    
end


for dd = 1:4
    if (dd ~= 1)
        path = 'Results/Empirical';
    else 
        path = 'Results/Simulation';
    end
    print_table_ESS(squeeze(mean_ESS_theta(:,:,:,dd)),...
         squeeze(mean_ESS_H(:,:,dd)),...
         squeeze(std_ESS_H(:,:,dd)),...
         squeeze(mean_AR_theta(:,:,dd)),...
         squeeze(mean_AR_H(:,dd)),...
         squeeze(mean_time(:,dd)),...
         params,DATA{dd},data_name{dd},path)
end



%% Times HMM adapt Different No of bins
N_bins = [10,20,30];
params = {'\mu','\phi','\sigma^2'};
DATA = {'', '_GSPC' , '_MSFT', '_IBM'};

Time_HMM_adapt = zeros(3,4);
for dd = 1:4
    data = DATA{dd};
    if (dd == 1)
        path = 'Results/Simulation/';
    else
        path = 'Results/Empirical/';
    end
    for jj = 1:3
        N_bin = N_bins(jj);
        clear -regexp ^time 
      
        load([path,'SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),data,'.mat'],...
            '-regexp','^time')
        aaa = who('-regexp','^time_');
        Time_HMM_adapt(jj,dd) = eval(char(aaa));
    end
end