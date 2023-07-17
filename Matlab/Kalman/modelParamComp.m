clear all; close all; fclose('all'); rng(0);
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

%% Parameters
% Monte Carlo parameters
simFactor=1;
N=3*24*simFactor; % points in time grid
M=10000000; % number of simulations

% Parameters for fish farming, see Table 4 Ewald2017
T=3; % time horizon
m=.1; % mortality rate
cr=1.1; % conversion rate
n0=10000; % number of recruits
% hc=3; % variable harvesting cost per kg (NOK/kg)
% fc=10; % variable feeding cost per kg per year (NOK/kg)
wInf=6; % asymptotic weight (kg)

% Bertalanffyo's growth function, see Footnote 20 Ewald2017
a=1.113;
b=1.097;
c=1.43;

r=0.0303;

gamma=0.0;

% Salmon
% mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
salmonParam1=[0.12, 0.23, 0.75, 2.6, 0.02, 0.01, 0.9, 0.57, 95];%down,down
salmonParam2=[0.12, 0.23, 0.75, 2.6, 0.02, 0.2, 0.9, 0.57, 95];%down,up
salmonParam3=[0.12, 0.23, 0.75, 2.6, 0.02, 0.6, 0.9, 0.57, 95];%up,up
salmonParam=salmonParam1;

hc=salmonParam(end).*.5.*.1;
fc=salmonParam(end).*.5.*.25;
salmonParam(end)=salmonParam(end).*.5+fc+hc;

% disp('Salmon parameters: mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0')
% disp(join(split(num2str(salmonParam),whitespacePattern),','))

% Soy
% mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0

soyParam1=[0.15, 0.5, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500];%low vol
soyParam2=[0.15, 1, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500];%medium vol
soyParam3=[0.15, 2, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500];%high vol
soyParam=soyParam1;

% disp('Soy parameters: mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0')
% disp(join(split(num2str(soyParam),whitespacePattern),','))

%% Simulate processes
%% Salmon
% Brownian motions, directly under Q
[W1,W2,t]=brownianMotions(T,N,M,salmonParam(7));

[salmonP,~]=schwartzTwoFactor(salmonParam1(8),salmonParam1(9),r,salmonParam1(2),salmonParam1(3),salmonParam1(4),salmonParam1(5),salmonParam1(6),W1,W2,t);
relP1=mean(salmonP./salmonParam1(9),2);
[salmonP,~]=schwartzTwoFactor(salmonParam2(8),salmonParam2(9),r,salmonParam2(2),salmonParam2(3),salmonParam2(4),salmonParam2(5),salmonParam2(6),W1,W2,t);
relP2=mean(salmonP./salmonParam1(9),2);
[salmonP,~]=schwartzTwoFactor(salmonParam3(8),salmonParam3(9),r,salmonParam3(2),salmonParam3(3),salmonParam3(4),salmonParam3(5),salmonParam3(6),W1,W2,t);
relP3=mean(salmonP./salmonParam1(9),2);

fig=newFigure();hold on;
plot(t,relP1,'r-');
plot(t,relP2,'g-');
plot(t,relP3,'b-');
legend('down,down','down,up','up,up','Orientation','horizontal','Location','southoutside')
title('Salmon Relative Mean Price')
exportgraphics(fig,'salmonRelMean.pdf')
clear W1 W2 salmonP;
%% Soy
% Brownian motions, directly under Q
ticBM=tic;
[W1,W2,t]=brownianMotions(T,N,M,soyParam(7));
ctimeBM=toc(ticBM);

% Schwartz 2-factor model under Q
[soyP,~]=schwartzTwoFactor(soyParam1(8),soyParam1(9),r,soyParam1(2),soyParam1(3),soyParam1(4),soyParam1(5),soyParam1(6),W1,W2,t);
relP1=mean(soyP./soyParam1(9),2);
[soyP,~]=schwartzTwoFactor(soyParam2(8),soyParam2(9),r,soyParam2(2),soyParam2(3),soyParam2(4),soyParam2(5),soyParam2(6),W1,W2,t);
relP2=mean(soyP./soyParam2(9),2);
[soyP,~]=schwartzTwoFactor(soyParam3(8),soyParam3(9),r,soyParam3(2),soyParam3(3),soyParam3(4),soyParam3(5),soyParam3(6),W1,W2,t);
relP3=mean(soyP./soyParam2(9),2);
% [soyP,soyDelta]=schwartzTwoFactor(soyParam(8),soyParam(9),r,soyParam(2),soyParam(3),soyParam(4),soyParam(5),soyParam(6),W1,W2,t);
% soyP=soyP./soyParam(9);

fig=newFigure();hold on;
plot(t,relP1,'r-');
plot(t,relP2,'g-');
plot(t,relP3,'b-');
legend('low vol','medium vol','high vol','Orientation','horizontal','Location','southoutside')
title('Soy Relative Mean Price')
exportgraphics(fig,'soyRelMean.pdf')
clear W1 W2 soyP;
