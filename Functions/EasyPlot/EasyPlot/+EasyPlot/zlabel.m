function zlabel(all_axes, label, varargin)
    if ~iscell(all_axes)
        zlabel(all_axes, label, varargin{:});
        return
    end

    for k = 1:size(all_axes, 1)
        for j = 1:size(all_axes, 2)
            zlabel(all_axes{k,j}, label, varargin{:});
        end
    end
end