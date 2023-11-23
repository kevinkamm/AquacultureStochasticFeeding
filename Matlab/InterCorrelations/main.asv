clear all; close all; fclose('all'); seed=0; rng(seed);
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
warning('off') %LSMC ill-conditioned matrices
%% Figure options
backgroundColor='w';
textColor='k';

%% Extract data from spreadsheet
r=0.0303;
T=3; % time horizon
frq=2; % every frq week/month
disp('Load data')

startYear=2018;
dataStartDate=datetime(startYear,1,1);
dataEndDate=datetime(startYear+T,1,1);

%% Select Panel
dataPanel=[[0,1,2,3,4,5]',[1,2,3,4,5,6]']./12; % small panel
% dataPanel=[[0,1,5,6,10,11]',1+[0,1,5,6,10,11]']./12; % mixed panel
% dataPanel=[[5,6,10,11,14,15]',1+[5,6,10,11,14,15]']./12; % medium panel

%% Joint Kalman
[paramsJoint,ss_attJoint,negLogLikelihoodJoint]=kalmanDaily2C(dataStartDate,dataEndDate,dataPanel,r,verbose=1);
[RhoJoint1,xRhoJoint1,LagJoint1]=corrFromFilter(ss_attJoint);
[RhoJoint2,xRhoJoint2,LagJoint2]=corrFromFilter(ss_attJoint,'differences');

%% Independent Kalman
[paramsSalmon,ss_attSalmon,negLogLikelihoodSalmon,datesSalmon]=kalmanDaily1C('Salmon',dataStartDate,dataEndDate,dataPanel,r,verbose=1);
[paramsSoy,ss_attSoy,negLogLikelihoodSoy,datesSoy]=kalmanDaily1C('Soy',dataStartDate,dataEndDate,dataPanel,r,verbose=1);
[C,iSoy,iSalmon]=intersect(datesSoy',datesSalmon');
ss_att=[ss_attSalmon(iSalmon,:),ss_attSoy(iSoy,:)];
[Rho1,xRho1,Lag1]=corrFromFilter(ss_att);
[Rho2,xRho2,Lag2]=corrFromFilter(ss_att,'differences');

%% Independent Cortazar
[paramsSalmonC,atSalmon,errSalmon,datesSalmon]=cortazarDaily1C('Salmon',dataStartDate,dataEndDate,r,verbose=1);
[paramsSoyC,atSoy,errSoy,datesSoy]=cortazarDaily1C('Soy',dataStartDate,dataEndDate,r,verbose=1);
[C,iSoy,iSalmon]=intersect(datesSoy',datesSalmon');
at=[atSalmon(iSalmon,:),atSoy(iSoy,:)];
[RhoC1,xRhoC1,LagC1]=corrFromFilter(at);
[RhoC2,xRhoC2,LagC2]=corrFromFilter(at,'differences');
