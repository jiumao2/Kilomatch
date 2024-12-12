function xlabel(all_axes, label)
    if ~iscell(all_axes)
        xlabel(all_axes, label);
        return
    end

    for k = 1:size(all_axes, 1)
        for j = 1:size(all_axes, 2)
            xlabel(all_axes{k,j}, label);
        end
    end
end