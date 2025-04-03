function setYLabelRow(all_axes, label, varargin)
    if ~iscell(all_axes)
        ylabel(all_axes, label, varargin{:});
        return
    end
    
    idx_leftmost = EasyPlot.Utils.getSideAxes(all_axes, 'left');
    ylabel(all_axes{idx_leftmost}, label, varargin{:});
end