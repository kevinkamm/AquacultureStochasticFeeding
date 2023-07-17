clear all; close all;
dataName='Soy';
data=load(['Daily/',dataName]);
prices=data.output;
dates=data.udates;
ttm=data.uTTM';

% filledPrices=fillmissing(prices,'nearest',2);

dataStartDate=datetime(2018,1,1);
dataEndDate=datetime(2022,12,31);
dateInd = dates >=dataStartDate & dates <=dataEndDate;
Fmarket=prices(dateInd,:);
days=dates(dateInd);
% yData=log(prices(dateInd,:));
% yData=log(filledPrices(dateInd,:));
r=0.0303;

%% Model Future Formula Cortazar-Schwartz Soy
paramsC=[r, 2.1975    0.45943    0.43818    2.2825    1.312     0.48937];
spot=load('savedSpots/spotCortazar');
aC=spot.aFin(2:end,:);

sigma11=paramsC(2);
sigma12=paramsC(3);
kappa1=paramsC(4);
alpha1=paramsC(5);
lambda1=paramsC(6);
rho12=paramsC(7);

tau1=ttm;

Atau1=(r-alpha1+lambda1./kappa1 + sigma12.^2./(2.*kappa1.^2)-sigma11.*sigma12.*rho12./kappa1).*tau1+...
sigma12.^2.*((1-exp(-2.*kappa1.*tau1))./(kappa1.^3) )./4+...
(alpha1.*kappa1-lambda1+sigma11.*sigma12.*rho12-sigma12.^2./kappa1).*((1-exp(-kappa1.*tau1))./(kappa1.^2));
Fcortazar=aC(:,1).*exp(-aC(:,2).*(1-exp(-kappa1.*tau1))./kappa1 + Atau1);
%% Model Future Formula Kalman Soy
% % Kalman, panel=[[1,2,3,4,5]',[2,3,4,5,6]']./12;
% paramsK=[0.1561921, 0.2128805, 0.2145218, 0.2967304, 0.02865752, 0.05314332, 0.5485235];
% spot=load('savedSpots/spotKalman1');

% Kalman, panel=[[5,10,15]',[6,11,16]']./12;
% paramsK=[0.1892143, 0.2601167, 0.2239957, 0.05631451, 0.9998518, 0.01000455, 0.8531076];
% spot=load('savedSpots/spotKalman2');

% Kalman, panel=[[0,1,5,6,10,11]',1+[0,1,5,6,10,11]']./12;
% paramsK=[0.1875447, 0.2110906, 0.1628161, 0.2174437, 0.2969554, 0.01000004, 0.7539584];
% spot=load('savedSpots/spotKalman3');

% Kalman, panel=[[0:20]',1+[0:20]']./12;
paramsK=[0.1944271, 0.2044593, 0.1275143, 0.1541256, 0.3355088, 0.01010092, 0.9065442];
spot=load('savedSpots/spotKalman4');



aK=spot.ss_att;aK(:,1)=exp(aK(:,1));
sigma11=paramsK(2);
sigma12=paramsK(3);
kappa1=paramsK(4);
alpha1=paramsK(5);
lambda1=paramsK(6);
rho12=paramsK(7);

tau1=ttm;

Atau1=(r-alpha1+lambda1./kappa1 + sigma12.^2./(2.*kappa1.^2)-sigma11.*sigma12.*rho12./kappa1).*tau1+...
sigma12.^2.*((1-exp(-2.*kappa1.*tau1))./(kappa1.^3) )./4+...
(alpha1.*kappa1-lambda1+sigma11.*sigma12.*rho12-sigma12.^2./kappa1).*((1-exp(-kappa1.*tau1))./(kappa1.^2));
Fkalman=aK(:,1).*exp(-aK(:,2).*(1-exp(-kappa1.*tau1))./kappa1 + Atau1);

%% Model Hybrid Soy
paramsH=(paramsC+paramsK)./2;
aH=(aK+aC)./2;
sigma11=paramsH(2);
sigma12=paramsH(3);
kappa1=paramsH(4);
alpha1=paramsH(5);
lambda1=paramsH(6);
rho12=paramsH(7);

tau1=ttm;

Atau1=(r-alpha1+lambda1./kappa1 + sigma12.^2./(2.*kappa1.^2)-sigma11.*sigma12.*rho12./kappa1).*tau1+...
sigma12.^2.*((1-exp(-2.*kappa1.*tau1))./(kappa1.^3) )./4+...
(alpha1.*kappa1-lambda1+sigma11.*sigma12.*rho12-sigma12.^2./kappa1).*((1-exp(-kappa1.*tau1))./(kappa1.^2));
Fhybrid=aK(:,1).*exp(-aK(:,2).*(1-exp(-kappa1.*tau1))./kappa1 + Atau1);

%%
dInd=size(Fkalman,1)-10;
fig=newFigure();
tl=tiledlayout(3,1);
%% Futures for one date
ax=nexttile; hold on;
plot(ttm,Fcortazar(dInd,:),'m--')
plot(ttm,Fkalman(dInd,:),'r--')
% plot(ttm,Fhybrid(dInd,:),'y--')
plot(ttm,Fmarket(dInd,:),'bx')
title(['Future prices for ',char(days(dInd+1))])
xlabel('Time to Maturity')
% legend(ax,'Cortazar','Kalman','Hybrid','Market','Location','southoutside','Orientation','horizontal')
legend(ax,'Cortazar','Kalman','Market','Location','southoutside','Orientation','horizontal')

%% Spots for all dates
ax=nexttile; hold on;
plot(aC(:,1),'mo')
plot(aK(:,1),'rx')
% plot(aH(:,1),'y.')
title('Spot for all dates')
xlabel('Dates')
% legend(ax,'Cortazar','Kalman','Hybrid','Location','southoutside','Orientation','horizontal')
legend(ax,'Cortazar','Kalman','Location','southoutside','Orientation','horizontal')

%% Convenience yield for all dates
ax=nexttile; hold on;
plot(aC(:,2),'mo')
plot(aK(:,2),'rx')
% plot(aH(:,2),'y.')
title('Convenience yield for all dates')
xlabel('Dates')
% legend(ax,'Cortazar','Kalman','Hybrid','Location','southoutside','Orientation','horizontal')
legend(ax,'Cortazar','Kalman','Location','southoutside','Orientation','horizontal')

exportgraphics(fig,'Figures/KalmanVsCortazar.pdf')
%%
fig=figure();hold on;
[X,Y]=meshgrid([1:size(Fkalman,1)],ttm);
% surf(X,Y,Fcortazar','EdgeColor','none','LineStyle','none','FaceColor','m'); hold on;
% surf(X,Y,Fkalman','EdgeColor','none','LineStyle','none','FaceColor','r'); hold on;
surf(X,Y,abs(Fkalman'-Fcortazar')./Fcortazar','EdgeColor','none','LineStyle','none'); hold on;
colorbar;
view(3);
xlabel('Dates')
ylabel('TTM')
zlabel('Price difference in %')
% legend('Cortazar','Kalman')