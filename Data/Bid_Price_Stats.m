clear all;
path = 'C:/Users/ancal/Desktop/DNB_Bid_Data_cleaned/';

% Ticks= {'AA', 'F','IBM' , 'JPM', 'KO', 'XRX'};
Ticks= {'IBM' , 'JPM', 'KO'};
Periods1 = {'2008103_9','20081010_10'};
X1 = zeros(9,6);
Periods2 = {'20100423_29','20100430_30'};
X2 = zeros(9,6);

for jj = 1:3
    tick = Ticks{jj};
    for ii = 1:2
        period = Periods1{ii};
        FileName = [path,'Data_',period,'_',tick,'_Bid.csv'];
        TickRet = xlsread(FileName,'B:B');
        Price = xlsread(FileName,'C:C');
        X1(1,2*(jj-1) + ii) = length(TickRet);
        X1(2,2*(jj-1) + ii) = mean(Price);
        X1(3,2*(jj-1) + ii) = mean(TickRet);
        X1(4,2*(jj-1) + ii) = std(TickRet);
        X1(5,2*(jj-1) + ii) = min(TickRet);
        X1(6,2*(jj-1) + ii) = max(TickRet);
        X1(7,2*(jj-1) + ii) = 100*sum(abs(TickRet)==0)/ length(TickRet);
        X1(8,2*(jj-1) + ii) = 100*sum(abs(TickRet)==1)/ length(TickRet);
        X1(9,2*(jj-1) + ii) = 100*sum((abs(TickRet)>=2) & (abs(TickRet)<=10))/length(TickRet);    
    end
end



for jj = 1:3
    tick = Ticks{jj};
    for ii = 1:2
        period = Periods2{ii};
        FileName = [path,'Data_',period,'_',tick,'_Bid.csv'];
        TickRet = xlsread(FileName,'B:B');
        Price = xlsread(FileName,'C:C');
        X2(1,2*(jj-1) + ii) = length(TickRet);
        X2(2,2*(jj-1) + ii) = mean(Price);
        X2(3,2*(jj-1) + ii) = mean(TickRet);
        X2(4,2*(jj-1) + ii) = std(TickRet);
        X2(5,2*(jj-1) + ii) = min(TickRet);
        X2(6,2*(jj-1) + ii) = max(TickRet);
        X2(7,2*(jj-1) + ii) = 100*sum(abs(TickRet)==0)/ length(TickRet);
        X2(8,2*(jj-1) + ii) = 100*sum(abs(TickRet)==1)/ length(TickRet);
        X2(9,2*(jj-1) + ii) = 100*sum((abs(TickRet)>=2) & (abs(TickRet)<=10))/length(TickRet);    
    end
end



% fname = '../../Paper/Tables/New/data_bid_short_comb.tex';
fname = 'data_bid_short_comb.tex';

FID = fopen(fname, 'w+');
% fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

% fprintf(FID, '\\begin{table} \n');
% fprintf(FID, '\\center \n');
fprintf(FID, '\\begin{singlespace} \n');
fprintf(FID, '\\begin{tabular}{lrrrrrr} \n');
fprintf(FID, '\\toprule \n');
 

fprintf(FID, ['\\multirow{2}{*}{October 2008} & \\multicolumn{2}{c}{IBM} ',...
    ' & \\multicolumn{2}{c}{KO}& \\multicolumn{2}{c}{JPM} \\\\ \n']);
fprintf(FID, ['\\cmidrule(r){2-3} \\cmidrule(r){4-5} \\cmidrule(r){6-7}  ',...
    ' & \\multicolumn{1}{c}{In} &  \\multicolumn{1}{c}{ Out} & \\multicolumn{1}{c}{In} ',...
    ' & \\multicolumn{1}{c}{Out} &  \\multicolumn{1}{c}{In} &  \\multicolumn{1}{c}{ Out} \\ \\midrule \n']);
fprintf(FID, ' Num. obs ' );
    for kk = 1:6
        fprintf(FID, '& %i', X1(1,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Avg. price ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X1(2,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Mean' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X1(3,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Std ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X1(4,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Min ' );
    for kk = 1:6
        fprintf(FID, '& %i', X1(5,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Max ' );
    for kk = 1:6
        fprintf(FID, '& %i', X1(6,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' \\%% 0  ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X1(7,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' \\%% $\\pm$ 1  ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X1(8,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' \\%% $\\pm$ 2--10  ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X1(9,kk));
    end
fprintf(FID,' \\\\ \\midrule \\midrule  \n');

 

fprintf(FID, ['\\multirow{2}{*}{April 2010} & \\multicolumn{2}{c}{IBM} ',...
    ' & \\multicolumn{2}{c}{KO}& \\multicolumn{2}{c}{JPM} \\\\ \n']);
fprintf(FID, ['\\cmidrule(r){2-3} \\cmidrule(r){4-5} \\cmidrule(r){6-7}  ',...
    ' & \\multicolumn{1}{c}{In} &  \\multicolumn{1}{c}{ Out} & \\multicolumn{1}{c}{In} ',...
    ' & \\multicolumn{1}{c}{Out} &  \\multicolumn{1}{c}{In} &  \\multicolumn{1}{c}{ Out} \\ \\midrule \n']);
fprintf(FID, ' Num. obs ' );
    for kk = 1:6
        fprintf(FID, '& %i', X2(1,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Avg. price ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X2(2,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Mean' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X2(3,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Std ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X2(4,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Min ' );
    for kk = 1:6
        fprintf(FID, '& %i', X2(5,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' Max ' );
    for kk = 1:6
        fprintf(FID, '& %i', X2(6,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' \\%% 0  ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X2(7,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' \\%% $\\pm$ 1  ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X2(8,kk));
    end
fprintf(FID,' \\\\  \n');

fprintf(FID, ' \\%% $\\pm$ 2--10  ' );
    for kk = 1:6
        fprintf(FID, '& %6.4f', X2(9,kk));
    end
fprintf(FID,' \\\\   \n');



fprintf(FID,'\\bottomrule \n');
fprintf(FID, '\\end{tabular}\n ');
fprintf(FID, '\\end{singlespace} \n');


% caption = ['\\caption{Posterior means, standard deviations (in parentheses) \n',...
%     'and the 2.5\\%%--97.5\\%% quantile ranges (in brackets) averaged over $50$ replications, \n',...
%     ' for $M=20,000$ posterior draws after a burn-in of $20,000$ \n' ,...
%     'for $T=20,000$ observations  generated from the Skellam and $\\Delta$NB models.}\n'];
% fprintf(FID, caption);

% label = ['\\label{tab:sim_res} \n'];
% fprintf(FID, label);

% fprintf(FID, '\\end{table}\n');

% fprintf(FID, '}');
fclose(FID);
