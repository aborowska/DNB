clear all;
% Ticks= {'AA', 'F','IBM' , 'JPM', 'KO', 'XRX'};
Ticks= {'IBM' , 'JPM', 'KO'};
Periods = {'2008','2010'};

for jj = 1:2
    for ii = 1:3
        period = Periods{jj};
        if strcmp(period,'2008')
            % ''' 2008 October '''
            InSampleStr = 'Data_2008103_9_';
            OutSampleStr = 'Data_20081010_10_';
            FullSampleStr = 'DataFull_2008103_10_';
        else
            % ''' 2010 April ''' 
            % DataFiles=['NYSE_201004_data_11.csv', 
            %            'NYSE_201004_data_12.csv', 
            %            'NYSE_201004_data_21.csv', 
            %            'NYSE_201004_data_22.csv', 
            %            'NYSE_201004_data_23.csv', 
            %            'NYSE_201004_data_24.csv', 
            %            'NYSE_201004_data_31.csv', 
            %            'NYSE_201004_data_32.csv', 
            %            'NYSE_201004_data_33.csv', 
            %            'NYSE_201004_data_34.csv', 
            %            'NYSE_201004_data_35.csv', 
            %            'NYSE_201004_data_36.csv', 
            %            'NYSE_201004_data_41.csv', 
            %            'NYSE_201004_data_51.csv', 
            %            'NYSE_201004_data_52.csv'] 
            % InSampleDays=[23, 26,27,28,29]
            % OutSampleDays=[30]

            InSampleStr = 'Data_20100423_29_';
            OutSampleStr = 'Data_20100430_30_';
            FullSampleStr = 'DataFull_20100423_30_';

        end

        %% DateTime
        filename = [FullSampleStr,Ticks{ii},'.csv'];
        deltaT_nonights = xlsread(filename,'H:H');

        if false
            [DateTime, Time] = ImportData(filename);

            t = datetime(DateTime,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            t1 = t(1:(end-1));
            t2 = t(2:end);
            t2(1) == t1(2)
            t2(1) == t(1)

            deltaT = milliseconds(t2-t1); % in milliseconds to maintain the precision
            deltaT = deltaT/1000; % in seconds
            T = length(deltaT); % 32432
            plot(deltaT)
            [B,I] = sort(deltaT);

            deltaT(I(end-4:end)); % 5 the largest - nights
            ind_break = sort(I(end-4:end));

            deltaT_nonights = deltaT;
            deltaT_nonights(ind_break) = [];
        end

        ff = figure(1);
            set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
            set(gcf,'defaulttextinterpreter','latex');
            plot(deltaT_nonights)
            hAxes = gca;
            hAxes.TickLabelInterpreter = 'latex';
            xlabel('Trades','FontSize',12)
            ylabel('Duration','FontSize',12)

            set(gcf,'PaperPositionMode','auto');
            plotname = ['duration_',period,'_',Ticks{ii},'.png'];
            print(ff,plotname,'-dpng','-r0')
            plotname = ['duration_',period,'_',Ticks{ii},'.eps'];
            print(ff,plotname,'-depsc','-r0')

        ff = figure(2);
            set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.5]);
            set(gcf,'defaulttextinterpreter','latex');
            hh = histogram(deltaT_nonights);
            BinEdges = hh.BinEdges;
            BinEdges = BinEdges(2:end);
            Values = hh.Values;
            HH = [BinEdges;Values]';
            bar(HH(:,1),log10(HH(:,2)),'FaceColor',[0    0.4470    0.7410])
            hAxes = gca;
            hAxes.TickLabelInterpreter = 'latex';
            xlabel('Duration','FontSize',12)
            ylabel('$\log_{10}$ frequency','FontSize',12)

            set(gcf,'PaperPositionMode','auto');
            plotname = ['frequency_',period,'_',Ticks{ii},'.png'];
            print(ff,plotname,'-dpng','-r0')
            plotname = ['frequency_',period,'_',Ticks{ii},'.eps'];
            print(ff,plotname,'-depsc','-r0')

        % save(['deltaT_hist_',Ticks{ii},'.mat'],'deltaT','deltaT_nonights','ind_break','HH');

    end
end



%%

for jj = 1:2
    for ii = 1:3
        period = Periods{jj};
        if strcmp(period,'2008')
            % ''' 2008 October '''
            InSampleStr = 'Data_2008103_9_';
            OutSampleStr = 'Data_20081010_10_';
            FullSampleStr = 'DataFull_2008103_10_';
        else
            % ''' 2010 April ''' 
            InSampleStr = 'Data_20100423_29_';
            OutSampleStr = 'Data_20100430_30_';
            FullSampleStr = 'DataFull_20100423_30_';

        end

        %% DateTime
   
        filename = [FullSampleStr,Ticks{ii},'.csv'];
        tick_return = xlsread(filename,'J:J'); 
        tick_return = sort(abs(tick_return));
        M = max(tick_return);
        counts = zeros(1+M,1);
        for tt = 0:M
            counts(tt+1) = sum(tick_return == tt);
        end
        counts = counts/length(tick_return);
        
        datafilename = ['counts_',period,'_',Ticks{ii},'.csv'];
        csvwrite(datafilename,counts);
    end
end



%% TradesIn
if false
    %        Calculating tick returns
    %         trades['TickReturn']=(trades['Price']/0.01).diff()   
    %         tradesIn=tradesIn[['TimeInSec', 'TickReturn','Price','LogReturn']]
    filename = [InSampleStr,Ticks{ii},'.csv'];
    DataIn = csvread(filename);
    Price = DataIn(:,3);
    TimeInSec = DataIn(:,1);
    TickReturn = DataIn(:,2);
    LogReturn = DataIn(:,4);

    PriceDiff = diff(Price/0.01);

    subplot(2,1,1)
    plot(PriceDiff)
    subplot(2,1,2)
    plot(TickReturn)
end