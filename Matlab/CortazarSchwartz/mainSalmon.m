clear all; close all; fclose('all'); rng(0);
try
    num_workers = str2num(getenv('SLURM_CPUS_PER_TASK'));
    old_threads = maxNumCompThreads(num_workers);
    pool=parpool('threads'); % multithreading
    fprintf('Chosen number of workers %d, Number of active workers %d\n',num_workers,pool.NumWorkers)
catch ME
end
pool=gcp('nocreate');
if isempty(pool)
%     pool=parpool('local'); % multiprocessing
    pool=parpool('threads'); % multithreading
end
availableGPUs = gpuDeviceCount('available');
if availableGPUs > 0
    gpuDevice([]); % clears GPU
    gpuDevice(1); % selects first GPU, change for multiple with spmd
end
%% Figure options
backgroundColor='w';
textColor='k';

%% Extract data from spreadsheet
disp('Load data')
dataName='Salmon';
data=load(['../Data/Daily/',dataName]);
prices=data.output;
dates=data.udates;
ttm=data.uTTM';

filledPrices=fillmissing(prices,'nearest',2);

dataStartDate=datetime(2018,1,1);
dataEndDate=datetime(2022,12,31);
dateInd = dates >=dataStartDate & dates <=dataEndDate;


yData=log(prices(dateInd,:));
% yData=log(filledPrices(dateInd,:));
r=0.0303;

%%
disp('Calibrate')
ticCal=tic;
[params,resOuter,aFin,resInner]=calibrate(r,ttm,yData);
ctimeCal=toc(ticCal);
fprintf('Calibration done in %g s with outer error %g and inner mean error %g\n',ctimeCal,resOuter,mean(resInner));

%%
[sigma11,sigma12,kappa1,alpha1,lambda1,rho,paramDesc]=paramUnpack(params);
paramTable=array2table(reshape(params,1,[]),'VariableNames',paramDesc);
disp(paramTable)
disp(join(split(num2str([nan,params(1:6),aFin(end,2),aFin(end,1)]),whitespacePattern),', '))
%%
fig=newFigure();
plot(filledPrices(dateInd,1),'bx')
plot(aFin(:,1),'r--')