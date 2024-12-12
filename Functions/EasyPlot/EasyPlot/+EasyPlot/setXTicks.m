function setXTicks(axes_all, ticks)
    ticks = sort(ticks);
    if ~iscell(axes_all)
        xticks(axes_all, ticks);
        return
    end

    for k = 1:size(axes_all,1)
        for j = 1:size(axes_all,2)
            xticks(axes_all{k,j}, ticks);
        end
    end
end