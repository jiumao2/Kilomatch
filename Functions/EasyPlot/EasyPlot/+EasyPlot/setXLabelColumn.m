function setXLabelColumn(all_axes, label, varargin)
    if ~iscell(all_axes)
        xlabel(all_axes, label, varargin{:});
        return
    end
    
    idx_bottom = EasyPlot.Utils.getSideAxes(all_axes, 'bottom');
    xlabel(all_axes{idx_bottom}, label, varargin{:});
end