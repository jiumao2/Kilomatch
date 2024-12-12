function setYLabelColumn(all_axes, label)
    if ~iscell(all_axes)
        ylabel(all_axes, label);
        return
    end
    
    for k = 1:length(all_axes)
        ylabel(all_axes{k},label);
    end
end