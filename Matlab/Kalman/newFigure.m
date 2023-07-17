function fig=newFigure(varargin)
backgroundColor='w';
textColor = 'k';
figureRatio = 'fullScreen';
visible='on';
if ~isempty(varargin)
for iV=1:2:length(varargin)
    switch varargin{iV}
        case 'backgroundColor'
            backgroundColor = varargin{iV+1};
        case 'figureRatio'
            figureRatio = varargin{iV+1};
        case 'textColor'
            textColor  = varargin{iV+1};
        case 'visible'
            visible  = varargin{iV+1};
    end
end
end
if strcmp(figureRatio,'square')
    fig=figure();
else
    fig=figure('units','normalized',...
               'outerposition',[0 0 1 1]);
end
hold on;
fig.Visible=visible;
if strcmp(fig.Visible,'on')
    % fig.WindowState = 'minimized';
end
fontsize=22;
linewidth=2;
markersize=12;
set(gca,'FontSize',fontsize)
set(gca,'defaultLineMarkerSize',markersize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
set(gca, 'color', backgroundColor);
set(gca, 'XColor', textColor);
set(gca, 'YColor', textColor);
set(gca, 'ZColor', textColor);
end
% dark blue: [0 0.4470 0.7410]
% orange: [0.8500 0.3250 0.0980]
% yellow: [0.9290 0.6940 0.1250]
% purple: [0.4940 0.1840 0.5560]
% green: [0.4660 0.6740 0.1880]
% light blue: [0.3010 0.7450 0.9330]
% red: [0.6350 0.0780 0.1840]