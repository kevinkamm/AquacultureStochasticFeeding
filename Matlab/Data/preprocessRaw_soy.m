clear all; close all;
%% Converts RAW data to Dayly data
% Output format:
% M=Maturity, TTM=Time-To-Maturity, P=Price
% Date / TTM &-0.1 & ... & 0.0      & 0.01      & ...
% Day 1      &     &     & P_1(0.0) & P_1(0.01) & ...
% Day 2      &     &     & P_2(0.0) & ....
% [...]
% If value for day does not exist its NaN
% For Soy negative TTM can happen
%
% Output will be splitted in 3 components: data matrix, date vector, TTM
% vector, such that data(i,j) corresponds to date(i) with TTM(j)

% Input
fileDir = 'Raw';
fileName = 'soybean futures 2003-2023';
fileExt = 'xlsx';
fileStr = [fileDir,'/',fileName,'.',fileExt];

% Output
saveDir = 'Daily'; mkDir([cd,'\',saveDir]);
saveName = 'Soy';
saveExt = 'mat';
saveStr = [saveDir,'/',saveName,'.',saveExt];

% Assumptions
delivery=32; % min(delivery,endOfMonth)
% delivery=15; % min(delivery,endOfMonth)
basis=0; %for year fraction, 0=actual
acc=1; % convert to 0.5 month precision for Fx naming


%% Process Raw Data
rawTable = readtable(fileStr,'VariableNamingRule','preserve');

%% Transform to same format as Salmon
stackedDates=[];
stackedDelivery=[];
stackedPrice=[];
for i=2:9:size(rawTable,2)
stackedDates=cat(1,stackedDates,rawTable{:,1});
stackedDelivery=cat(1,stackedDelivery,datetime(rawTable{:,i},'InputFormat','MMM yy','Locale','en_US'));
tmpBid=rawTable{:,i+1};
if iscell(tmpBid)
    tmpBid = [cellfun(@str2double,tmpBid)];
end
tmpAsk=rawTable{:,i+2};
if iscell(tmpAsk)
    tmpAsk = [cellfun(@str2double,tmpAsk)];
end
stackedPrice=cat(1,stackedPrice,(tmpBid+tmpAsk)./2);
end
stackedDates.Format='yyyy-MM-dd';
stackedDates.Format='yyyy-MM-dd';
[stackedDates,sorted]=sort(stackedDates,1);
stackedDelivery=stackedDelivery(sorted,:);
stackedPrice=stackedPrice(sorted,:);
stackedYears=year(stackedDelivery);
stackedMonths=month(stackedDelivery);
stackedTable=table(stackedDates,stackedYears,stackedMonths,stackedPrice);
ind=~isnat(stackedTable{:,1});
dates=stackedTable{ind,1};
rawYears=stackedTable{ind,2};
rawMonths=stackedTable{ind,3};
rawPrice=stackedTable{ind,4};
maturity=dates;

%% Use same code as for Salmon
rawPrice(rawPrice<=0)=nan;
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
function x=myround(x,acc)
    x=round(x./acc).*acc;
end