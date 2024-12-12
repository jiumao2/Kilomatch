function title(all_axes, text, varargin)
    if ~iscell(all_axes)
        title(all_axes, text, varargin{:});
        return
    end

    for k = 1:size(all_axes, 1)
        for j = 1:size(all_axes, 2)
            title(all_axes{k,j}, text, varargin{:});
        end
    end
end