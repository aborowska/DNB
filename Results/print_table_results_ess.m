function print_table_results_ess(param_true, Means_DNB, Std_DNB, Q025_DNB, Q975_DNB, ...
    Means_IF, Std_IF)

param = {'$\mu$','$\varphi$','$\sigma^2_{\eta}$','$\gamma$','$\beta_1$','$\beta_2$','$\nu$'};

    fname = ['DNB_sim_results_ess.tex'];
 
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cccc ccc} \n');
    fprintf(FID, '\\toprule \n');
 
     fprintf(FID, ['Parameter && True  & Mean & Std & 95\\%% range && Mean IF & Std IF \\\\ \\cline{1-1}  \\cline{3-6} \\cline{8-9}\n']);
   
    for ii = 1:7
        fprintf(FID, '%s  ', param{ii}); 
        fprintf(FID,' &&  %6.4f ', param_true(ii)); 
 
        fprintf(FID,' &  %6.4f ', mean(Means_DNB(ii,:))); 
 
        fprintf(FID,' &  (%6.4f) ', mean(Std_DNB(ii,:)));
   
        fprintf(FID,' &  [%6.4f,%6.4f] ', mean(Q025_DNB(ii,:)), mean(Q975_DNB(ii,:)));
 
        fprintf(FID,' &&  %6.4f ', Means_IF(ii)); 
 
        fprintf(FID,' &  (%6.4f) ', Std_IF(ii));
            
        fprintf(FID,' \\\\  \n');
    end
    
    fprintf(FID,'\\bottomrule \n');
 

%   fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{Average posterior standard deviations in parantheses '...
%         'at lag equal to the ' ...
%         'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{AR: acceptance rate ' ...
%     'of the MH RW algorithm. (for $h_{t}$ average over imputations)}}  \\\\ \n']);
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Posterior means, standard deviations (in parentheses), \n',...
        ' the 2.5\\%%--97.5\\%% quantile ranges (in brackets) and ',...
        ' the mean inefficiency factors (IF) and their standard deviations ',...
        ' averaged over MC $50$ replications, \n',...
        ' for $M=20,000$ posterior draws after a burn-in of $20,000$ \n' ,...
        'for $T=20,000$ observations  generated from the $\\Delta$NB models.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:sim_res} \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end