function idx = getSideAxes(all_axes, position)
    pos_left = zeros(1,length(all_axes));
    pos_height = zeros(1,length(all_axes));
    for k = 1:length(all_axes)
        pos_left(k) = all_axes{k}.Position(1);
        pos_height(k) = all_axes{k}.Position(2);
    end
    if strcmpi(position,'left') || strcmpi(position,'l')
        [~, idx] = min(pos_left);
    elseif strcmpi(position,'right') || strcmpi(position,'r')
        [~, idx] = max(pos_left);
    elseif strcmpi(position,'bottom') || strcmpi(position,'b')
        [~, idx] = min(pos_height);
    elseif strcmpi(position,'top') || strcmpi(position,'t')
        [~, idx] = max(pos_height);
    end
end