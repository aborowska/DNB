function fname = print_table_ESS(ESS_theta, ESS_H, std_ESS_H,...
    AR_theta, AR_H, time,...
    params,data,data_name,path)

    fname = [path,'/ESS',data,'.tex'];
    lags = {'sig','40', '100','1000'};
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{ccc cc cc} \n');
    fprintf(FID, '\\hline \n');
    fprintf(FID, [' & ESS lag&& DA KSC & DA RW & HMM fix & HMM adapt \\\\ \\hline  \\hline\n']);
    
    fprintf(FID, ' \\multicolumn{2}{c}{time (s)}&');
    for ii = 1:4
        fprintf(FID,' & %6.4f ', time(ii)); 
    end
    fprintf(FID,' \\\\  \\hline \n');

    fprintf(FID, ['\\multicolumn{7}{c}{$\\theta$} \\\\ \\hline \n']);  
    for pp = 1:3
        fprintf(FID,'\\multirow{4}{*}{$%s$}  ',params{pp});
        for jj = 1:4
            fprintf(FID,' & %s & ',lags{jj});
            for ii = 1:4
                fprintf(FID,' & %6.4f ', ESS_theta(jj,ii,pp)); 
            end
            fprintf(FID,' \\\\ \n');        
        end
        
        fprintf(FID,' & AR & & -- ');
        for ii = 2:4
            fprintf(FID,' & %6.4f ', AR_theta(ii,pp)); 
        end
        fprintf(FID,' \\\\ [1.3ex] \n');        
    end
    fprintf(FID, '\\hline \n');

    fprintf(FID, ['\\multicolumn{7}{c}{$ \\bm{h} $} \\\\ \\hline \n']);  
    for jj = 1:4
        fprintf(FID,'mean & %s & ',lags{jj});
        for ii = 1:4
            fprintf(FID,' & %6.4f ', ESS_H(jj,ii)); 
        end
        fprintf(FID,' \\\\ \n'); 
        
        fprintf(FID,'std & %s & ',lags{jj});
        for ii = 1:4
            fprintf(FID,' & [%6.4f] ', std_ESS_H(jj,ii)); 
        end
        fprintf(FID,' \\\\  [1ex]\n');           
    end     
    
    fprintf(FID,' & AR & ');
    for ii = 1:4
        fprintf(FID,' & %6.4f ', AR_H(ii)); 
    end
    fprintf(FID,' \\\\  [1ex]\n');           
        
    fprintf(FID, '\\hline \n');
 
    fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{sig: lag equal to the ' ...
    'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
    fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{AR: acceptance rate ' ...
    'of the MH RW algorithm.}}  \\\\ \n']);
    fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{mean/std: over the imputed states.}}  \\\\ \n']);
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Effective sample sizes (ESS) at different lags ',...
        ' for $M=10000$ posterior draws after a burn-in of $10000$ ' ,...
        'for $T=2000$ observations of \\textbf{',data_name ,'} data.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:ESS',data,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end