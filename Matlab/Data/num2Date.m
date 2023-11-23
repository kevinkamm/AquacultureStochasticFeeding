function num2Date(ax,dates)
    ticks=ax.XTickLabel;
    nDates=length(dates);
    for ti = 1:length(ticks)
        curr=ticks(ti);
        num=str2num(curr{1})+1;
        if num <= nDates
            ticks(ti)={string(dates(num))};
        else
            ticks(ti)={" "};
        end
    end
    ax.XTickLabel=ticks;
end