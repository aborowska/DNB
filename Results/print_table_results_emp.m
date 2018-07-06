function print_table_results_emp(Data,Mean_emp, Std_emp, Q025_emp, Q975_emp, IF_emp)
    
    Data_latex = strrep(Data,'_',' ');
    
    param = {'$\mu$','$\varphi$','$\sigma^2_{\eta}$','$\gamma$','$\beta_1$','$\beta_2$','$\nu$'};

    fname = [Data,'_results.tex'];
 
    
    FID = fopen(fname, 'w+');
%     fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\tabcolsep=0.07cm \n');
    fprintf(FID, '\\begin{footnotesize} \n');
    fprintf(FID, '\\begin{singlespace} \n');
    
    fprintf(FID, '\\begin{tabular}{ccc cccc } \n');
    fprintf(FID, '\\toprule \n');
 
     fprintf(FID, ['Parameter &&& OrdN  & Ord$t$ & Sk & $\\Delta$NB \\\\ \\hline\n']);
   
    for ii = 1:7
        if (ii < 7)
    %         fprintf(FID, '\\multirow{5}{*}{%s}   ', param{ii}); 
            fprintf(FID, ' %s  ', param{ii}); 

            fprintf(FID, ' & Mean   & ' );          
            for jj = 1:4
                fprintf(FID,' &  %6.4f ', Mean_emp(ii,jj)); 
            end
            fprintf(FID,' \\\\  \n');

            fprintf(FID, ' & Std   & ' );         
            for jj = 1:4       
                fprintf(FID,' &  (%6.4f) ', Std_emp(ii,jj));
            end
            fprintf(FID,' \\\\  \n');

            fprintf(FID, ' & 95\\%%  & ' );           
            for jj = 1:4
                fprintf(FID,' &  [%6.4f,%6.4f] ', Q025_emp(ii,jj), Q975_emp(ii,jj));
            end
            fprintf(FID,' \\\\  \n');

            fprintf(FID, ' & IF  & ' );                 
            for jj = 1:4
                fprintf(FID,' &  %6.4f ', IF_emp(ii,jj)); 
            end
            fprintf(FID,' \\\\ [1.0ex] \n');     
        else
            fprintf(FID, ' %s  ', param{ii}); 

            fprintf(FID, ' & Mean   & ' );          
            for jj = 2:2:4
                fprintf(FID,' &&  %6.4f ', Mean_emp(ii,jj)); 
            end
            fprintf(FID,' \\\\  \n');

            fprintf(FID, ' & Std   & ' );         
            for jj = 2:2:4
                fprintf(FID,'& &  (%6.4f) ', Std_emp(ii,jj));
            end
            fprintf(FID,' \\\\  \n');

            fprintf(FID, ' & 95\\%%  & ' );           
            for jj = 2:2:4
                fprintf(FID,' &&  [%6.4f,%6.4f] ', Q025_emp(ii,jj), Q975_emp(ii,jj));
            end
            fprintf(FID,' \\\\  \n');

            fprintf(FID, ' & IF  & ' );                 
            for jj = 2:2:4
                fprintf(FID,' &&  %6.4f ', IF_emp(ii,jj)); 
            end
            fprintf(FID,' \\\\ [1.0ex] \n');                 
        end
    end
    
    fprintf(FID,'\\bottomrule \n');
 

%   fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{Average posterior standard deviations in parantheses '...
%         'at lag equal to the ' ...
%         'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{AR: acceptance rate ' ...
%     'of the MH RW algorithm. (for $h_{t}$ average over imputations)}}  \\\\ \n']);
    fprintf(FID, '\\end{tabular}\n ');
    fprintf(FID, '\\end{singlespace}\n ');
    fprintf(FID, '\\end{footnotesize}\n ');
    caption = ['\\caption{Posterior means, standard deviations (Std, in parentheses), \n',...
        ' the 2.5\\%%--97.5\\%% quantile ranges (95\\%%, in brackets) and ',...
        ' the inefficiency factors (IF) ',...
        ' for $M=20,000$ posterior draws after a burn-in of $20,000$ \n' ,...
        'for ',Data_latex,' data.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:sim_res_',Data,'} \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
%     fprintf(FID, '}');
    fclose(FID);
end