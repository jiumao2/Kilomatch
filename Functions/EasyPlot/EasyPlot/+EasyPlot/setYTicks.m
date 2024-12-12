function setYTicks(axes_all, ticks)
    ticks = sort(ticks);
    if ~iscell(axes_all)
        yticks(axes_all, ticks);
        return
    end

    for k = 1:size(axes_all,1)
        for j = 1:size(axes_all,2)
            yticks(axes_all{k,j}, ticks);
        end
    end
end