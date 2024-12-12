function setZTicks(axes_all, ticks)
    ticks = sort(ticks);
    if ~iscell(axes_all)
        zticks(axes_all, ticks);
        return
    end

    for k = 1:size(axes_all,1)
        for j = 1:size(axes_all,2)
            zticks(axes_all{k,j}, ticks);
        end
    end
end