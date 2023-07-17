clear all; close all;
% comdty='Salmon';
comdty='Soy';
comdtyFull=load(['Daily\',comdty]);

dates=comdtyFull.udates;
ttm=comdtyFull.uTTM';
price=comdtyFull.output;

ind=ttm<=5;
[X,Y]=meshgrid([1:length(dates)]',ttm(ind));
filledPrice = fillmissing(price(:,ind),'nearest',2);
fig=newFigure();
plt=surf(X,Y,filledPrice','EdgeColor','none','LineStyle','none');
xlabel('date')
ylabel('time-to-maturity')
zlabel('price')
view(3);
exportgraphics(fig,['Figures/',comdty,'MarketSurface.pdf'])
% view([0,-1,0])


% fig=figure();hold on;
% plot(ttm(ind),filledPrice(500:502,ind)','-')
% fig=figure();hold on;
% plot(ttm(ind),filledPrice(1000:1002,ind)','-')
% fig=figure();hold on;
% plot(ttm(ind),filledPrice(1500:1502,ind)','-')
% fig=figure();hold on;
% plot(ttm(ind),filledPrice(2000:2002,ind)','-')
% fig=figure();hold on;
% plot(ttm(ind),filledPrice(2500:2502,ind)','-')