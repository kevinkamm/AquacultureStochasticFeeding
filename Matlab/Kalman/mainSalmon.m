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
ttm=data.uTTM;

dataStartDate=datetime(2018,1,1);
dataEndDate=datetime(2022,12,31);
dataPanel=[[0,1,2,3,4,5]',[1,2,3,4,5,6]']./12; % small panel
% dataPanel=[[0,1,5,6,10,11]',1+[0,1,5,6,10,11]']./12; % mixed panel
% dataPanel=[[5,6,10,11,14,15]',1+[5,6,10,11,14,15]']./12; % medium panel
% dataPanel=[[10,11,14,15,19,20]',1+[10,11,14,15,19,20]']./12; % large panel
dateInd = dates >=dataStartDate & dates <=dataEndDate;
ttmInd = ttm >= 0 & ttm <= dataPanel(end);

filledPrices=fillmissing(prices,'nearest',2);

yData=zeros(sum(dateInd),size(dataPanel,1));
panel=zeros(1,size(dataPanel,1));
for p=1:size(dataPanel,1)
    indP= ttm>dataPanel(p,1) & ttm<=dataPanel(p,2);
    yData(:,p)=mean(filledPrices(dateInd,indP),2);
    panel(p)=mean(ttm(indP));
end


yData=log(yData);
indContracts=1:size(yData,2); % start from 2 to exclude spot
r=0.0303;

% dt = 1/365; % daily time step
dt = years(mean(diff(dates(dateInd)))); % =250 daily time step
% dt=1/300;

%%
disp('Calibrate')
% varEpsMode='single';
varEpsMode='diag';
% varEpsMode='tridiag';
ticEM=tic;
[params,negLogLikelihood]=EMschwartz(yData(:,indContracts),r,dt,panel(indContracts),varEpsMode);
ctimeEM=toc(ticEM);
fprintf('Elapsed time for EM-Algo: %g s with log-likelihood %g\n',ctimeEM,-negLogLikelihood)

[P0,a0,c,T,varEta,dtau,Ztau,varEpsMode,~]=schwartzSSM(params,yData(:,indContracts),r,dt,panel(indContracts),varEpsMode);
[logL,ss_att]=kalmanFilter(P0,a0,c,T,varEta,dtau,Ztau,varEpsMode,yData(:,indContracts));

%%
paramDesc=cell(length(params),1);
paramDesc{1}='\mu';paramDesc{2}='\sigma_1';paramDesc{3}='\sigma_2';
paramDesc{4}='\kappa';paramDesc{5}='\alpha';paramDesc{6}='\lambda';
paramDesc{7}='\rho';
for i=1:length(paramDesc)-7
    paramDesc{7+i}=sprintf('s%d',i);
end
paramTable=array2table(params,'VariableNames',paramDesc);
disp(paramTable)

%%
fig=newFigure();

plot(exp(yData(2:end,1)),'k','linewidth',1);
plot(exp(yData(2:end,2:end)),'r--','linewidth',1);
plot(exp(ss_att(:,1)),'b','linewidth',1);

%%
% delta0=r is a good enough approx
% disp('Soy parameters: mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, Price0')
disp(join(split(num2str([params(1:7),ss_att(end,2),exp(ss_att(end,1))]),whitespacePattern),', '))

% %%
% ind1=1:length(panel);%only one at the moment
% a0(1)=exp(a0(1));
% [err,S10,delta10]=measureModelFit(a0,params,r,panel,yData,ind1);
% disp(mean(err))
% fig=newFigure();
% plot(S10,'bx');
% plot(exp(yData(:,1)),'k-');
% plot(exp(yData(:,2:end)),'r--');
% 
% %%
% fig=newFigure();
% subplot(2,1,1);hold on;
% plot(exp(ss_att(:,1)))
% plot(S10)
% 
% subplot(2,1,2);hold on;
% plot(ss_att(:,2))
% plot(delta10)