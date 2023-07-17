clear all; close all;
%% Converts RAW data to Dayly data
% Output format:
% M=Maturity, TTM=Time-To-Maturity, P=Price
% Date / TTM &-0.1 & ... & 0.0      & 0.01      & ...
% Day 1      &     &     & P_1(0.0) & P_1(0.01) & ...
% Day 2      &     &     & P_2(0.0) & ....
% [...]
% If value for day does not exist its NaN
% For Salmon negative TTM can happen
%
% Output will be splitted in 3 components: data matrix, date vector, TTM
% vector, such that data(i,j) corresponds to date(i) with TTM(j)

% Input
fileDir = 'Raw';
fileName = 'Forwardprices_full';
fileExt = 'csv';
fileStr = [fileDir,'/',fileName,'.',fileExt];

% Output
saveDir = 'Daily'; mkDir([cd,'\',saveDir]);
saveName = 'Salmon';
saveExt = 'mat';
saveStr = [saveDir,'/',saveName,'.',saveExt];

% Assumptions
% assume Fx means max x-month maturity
delivery=32; % min(delivery,endOfMonth)
basis=0; %for year fraction, 0=actual
acc=1; % convert to 0.5 month precision for Fx naming


%% Process Raw Data
rawTable = readtable(fileStr,'VariableNamingRule','preserve','Delimiter',';','DecimalSeparator',',');

dates=rawTable{:,1};
rawYears=rawTable{:,2};
rawMonths=rawTable{:,3};
rawPrice=rawTable{:,4};
rawPrice(rawPrice<=0)=nan;
maturity=dates;

maturity.Year=rawYears;
maturity.Month=rawMonths;
maturity.Day=deliveryRule(delivery,rawMonths,rawYears);
ttm=yearfrac(dates,maturity,basis);

startDate = dates(1);
endDate = dates(end);

% fdates = startDate : endDate; %full days
% bdates = busdays(startDate,endDate);% business days
udates = unique(dates); %unique given days

uTTM = unique(ttm);

output = nan.*ones(length(udates),length(uTTM));
% for idx = 1:numel(output)
%     output{idx} = nan;
% end

for d=1:length(udates)
    ind=dates==udates(d);
    currDates=dates(ind);
    currMat=maturity(ind);
    currTTM=ttm(ind);
    currPrice=rawPrice(ind);

    indTTM=logical(sum(uTTM==currTTM',2));

    output(d,indTTM)=currPrice;
end
save(saveStr,"udates");
save(saveStr,"uTTM","-append");
save(saveStr,"output","-append");


%% Auxiliary
function flag=mkDir(dir)
    if ~exist(dir,"dir")
        flag=mkdir(dir);
    else
        flag=1;
    end
end
function day=deliveryRule(delivery,month,year)
    day=min(delivery,eomday(year,month));
end