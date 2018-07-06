function print_table_results(param_true, ChainsDNB, ChainsSk)

param = {'$\mu$','$\varphi$','$\sigma^2_{\eta}$','$\gamma$','$\beta_1$','$\beta_2$','$\nu$'};

Means_DNB = squeeze(mean(ChainsDNB,1));
Std_DNB = squeeze(std(ChainsDNB,1));
Q025_DNB = squeeze(quantile(ChainsDNB,0.025,1));
Q975_DNB = squeeze(quantile(ChainsDNB,0.975,1));
 
Means_Sk = squeeze(mean(ChainsSk,1));
Std_Sk = squeeze(std(ChainsSk,1));
Q025_Sk = squeeze(quantile(ChainsSk,0.025,1));
Q975_Sk = squeeze(quantile(ChainsSk,0.975,1));
 

     
    fname = ['DNB_sim_results.tex'];
%     lags = {'sig','40', '100','1000'};
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{c cc cc cc} \n');
    fprintf(FID, '\\toprule \n');
%     fprintf(FID, [' & ESS lag&& DA KSC & DA RW & HMM fix & HMM adapt \\\\ \\hline  \\hline\n']);
    fprintf(FID, ['Parameter && True  && Skellam && $\\Delta$NB \\\\ \\hline  \\hline\n']);
    
    for ii = 1:7
        fprintf(FID, '%s  ', param{ii}); 
        fprintf(FID,' &&  %6.4f ', param_true(ii)); 
        if ii == 7
            fprintf(FID,' &&   '); 
        else
            fprintf(FID,' &&  %6.4f ', mean(Means_Sk(ii,:))); 
        end
        fprintf(FID,' &&  %6.4f ', mean(Means_DNB(ii,:)));
        fprintf(FID,' \\\\  \n');

        fprintf(FID,' &&  '); 
        if ii == 7
            fprintf(FID,' &&   ');
        else
            fprintf(FID,' &&  (%6.4f) ', mean(Std_Sk(ii,:)));
        end
        fprintf(FID,' &&  (%6.4f) ', mean(Std_DNB(ii,:)));
        fprintf(FID,' \\\\  \n');
        
        fprintf(FID,' &&  '); 
        if ii == 7
            fprintf(FID,' &&  ');
        else
            fprintf(FID,' &&  [%6.4f,%6.4f] ', mean(Q025_Sk(ii,:)), mean(Q975_Sk(ii,:)));
        end
        fprintf(FID,' &&  [%6.4f,%6.4f] ', mean(Q025_DNB(ii,:)), mean(Q975_DNB(ii,:)));
        fprintf(FID,' \\\\ [1ex] \n');
    end
    fprintf(FID,'\\bottomrule \n');
 

%   fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{Average posterior standard deviations in parantheses '...
%         'at lag equal to the ' ...
%         'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{AR: acceptance rate ' ...
%     'of the MH RW algorithm. (for $h_{t}$ average over imputations)}}  \\\\ \n']);
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Posterior means, standard deviations (in parentheses) \n',...
        'and the 2.5\\%%--97.5\\%% quantile ranges (in brackets) averaged over $50$ replications, \n',...
        ' for $M=20,000$ posterior draws after a burn-in of $20,000$ \n' ,...
        'for $T=20,000$ observations  generated from the Skellam and $\\Delta$NB models.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:sim_res} \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end