%%
filename = 'H:\Desktop\DNBFullData\Data_2008103_9_IBM.csv'; 
DataIBM_In = csvread(filename);
N_In = size(DataIBM_In,1);
TickReturn_In = DataIBM_In(:,2);
Price_In = DataIBM_In(:,3);

plot(Price_In - Price(2:N_In+1));
TickReturn_comp = round(diff(Price/0.01));
TickReturn_In_comp = [0; round(diff(Price_In/0.01))];


TickReturn_skip = TickReturn(2:end);
DiffTickRet = TickReturn_skip - TickReturn_comp;
ind = find(DiffTickRet~=0); 
DiffTickRet(ind)
%%
filename = 'H:\Desktop\DNBFullData\DataFull_2008103_10_IBM.csv';
delimiter = ',';
formatSpec = '%*s%*s%*s%s%s%s%*s%*s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
Price = cell2mat(raw(:, 1));
BidPrice = cell2mat(raw(:, 2));
AskPrice = cell2mat(raw(:, 3));
TimeInSec = cell2mat(raw(:, 4));
TickReturn = cell2mat(raw(:, 5));

%%
TickReturn  = diff(Price/0.01);
TickReturnAsk  = [0; diff(AskPrice(2:N_In+1)/0.01)];
TickReturnBid  = [0; diff(BidPrice(2:N_In+1)/0.01)];

TickReturn = TickReturn(1:68002);
TickReturnAsk = TickReturnAsk(1:68002);
TickReturnBid = TickReturnBid(1:68002);