function setXLabelColumn(all_axes, label)
    if ~iscell(all_axes)
        xlabel(all_axes, label);
        return
    end
    
    idx_bottom = EasyPlot.Utils.getSideAxes(all_axes,'bottom');
    xlabel(all_axes{idx_bottom},label);
end