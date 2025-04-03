function setYLabelColumn(all_axes, label, varargin)
    if ~iscell(all_axes)
        ylabel(all_axes, label, varargin{:});
        return
    end
    
    for k = 1:length(all_axes)
        ylabel(all_axes{k},label, varargin{:});
    end
end