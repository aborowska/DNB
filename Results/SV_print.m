clear all
EXT = {'DA','DA_RW_eff','HMM_eff_Nbin30','HMM_adapt_eff_Nbin30'};
METH = {'DA KSC', 'DA RW', 'HMM fix','HMM adapt'};
params = {'\mu','\phi','\sigma^2'};
DATA = {'', '_GSPC' , '_MSFT', '_IBM'};
data_name = {'simulated','S\\&P500','MSFT','IBM'};


TH_mean = zeros(4,3,4);
TH_std = zeros(4,3,4);

ind_sel = [200,600,1000,1400,1800];
ind_sel2 = [1,5,9,13,17];

HH_mean_sel = zeros(4,5,4);
HH_std_sel = zeros(4,5,4);

TH_ESS = zeros(4,3,4);
HH_ESS_sel = zeros(4,5,4);


TH_ESS = zeros(4,3,4);
HH_ESS_sel = zeros(4,5,4);


TH_AR = zeros(4,3,4);
HH_AR_mean = zeros(4,4);

Time = zeros(4,4);

for dd = 1:4 % DATA
    data = DATA{dd};
    if ~isempty(char(data))
        path_results = 'Results/Empirical';
    else 
        path_results = 'Results/Simulation';
    end

    for jj = 1:4 % METHOD
        ext = EXT{jj};
        clear -regexp ^ESS_ ^time ^A_ ^theta ^H_ 
        load([path_results,'/SV_results_param_',ext,char(data),'.mat'],...
            '-regexp','^ESS_\w*sig','^time','^A_','^theta_','^H_subset')
        clear theta_true
        aaa = who('-regexp','^ESS_theta_');
        bbb = who('-regexp','^ESS_H_');
        ccc = who('-regexp','^A_theta_');
        ddd = who('-regexp','^A_H_');        
        eee = who('-regexp','^time_'); 
        fff = who('-regexp','^theta_');
        ggg = who('-regexp','^H_subset_');
                
        TH_ESS(jj,:,dd) = eval(char(aaa));
        xxx = eval(char(bbb));
        if (jj < 3)
            HH_ESS_sel(jj,:,dd) = xxx(ind_sel);
        else
            HH_ESS_sel(jj,:,dd) = xxx(ind_sel/2);                
        end
        TH_mean(jj,:,dd) = mean(eval(char(fff)));
        TH_std(jj,:,dd) = std(eval(char(fff)));

        xxx = eval(char(ggg));
        HH_mean_sel(jj,:,dd) = mean(xxx(:,ind_sel2));
        HH_std_sel(jj,:,dd) = std(xxx(:,ind_sel2)); 
        
        if (jj > 1)
            TH_AR(jj,:,dd) = mean(eval(char(ccc)));
        end
        HH_AR_mean(jj,dd) = mean(eval(char(ddd)));
        Time(jj,dd) = eval(char(eee));       
    end  
    
end


for dd = 1:4
    if (dd ~= 1)
        path = 'Results/Empirical';
    else 
        path = 'Results/Simulation';
    end
    print_table_results(squeeze(TH_mean(:,:,dd)),...
         squeeze(TH_std(:,:,dd)),...
         squeeze(TH_ESS(:,:,dd)),...
         squeeze(TH_AR(:,:,dd)),...
         squeeze(HH_mean_sel(:,:,dd)),...
         squeeze(HH_std_sel(:,:,dd)),...
         squeeze(HH_ESS_sel(:,:,dd)),...
         squeeze(HH_AR_mean(:,dd)),...
         Time(:,dd),...
         DATA{dd},data_name{dd},path)
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