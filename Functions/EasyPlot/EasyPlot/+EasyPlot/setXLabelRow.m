function setXLabelRow(all_axes, label)
    if ~iscell(all_axes)
        xlabel(all_axes, label);
        return
    end

    for k = 1:length(all_axes)
        xlabel(all_axes{k},label);
    end
end