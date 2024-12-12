function setYLabelRow(all_axes, label)
    if ~iscell(all_axes)
        ylabel(all_axes, label);
        return
    end
    
    idx_leftmost = EasyPlot.Utils.getSideAxes(all_axes,'left');
    ylabel(all_axes{idx_leftmost},label);
end