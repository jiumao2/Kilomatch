function xLimOut = setXLim(all_axes,mode)
    if nargin<2
        mode = 'largest';
    end
    xLimOut = nan(1,2);

    if ~iscell(all_axes)
        if ~isnumeric(mode) || length(mode)~=2
            error('Wrong mode!')
        else
            xLimOut = mode;
        end
        set(all_axes,'XLim',xLimOut);
        return
    end

    if isnumeric(mode)
        if length(mode) == 1
            xLimOut = all_axes{mode}.XLim;
        elseif length(mode) == 2
            xLimOut = mode;
        else
            error('Wrong mode!');
        end
    elseif ischar(mode) || isstring(mode)
        if strcmpi(mode,'largest')
            xLimLeft = Inf;
            xLimRight = -Inf;
            for k = 1:size(all_axes,1)
                for j = 1:size(all_axes,2)
                    left = all_axes{k,j}.XLim(1);
                    right = all_axes{k,j}.XLim(2);
                    xLimLeft = min(xLimLeft,left);
                    xLimRight = max(xLimRight,right);
                end
            end
            xLimOut = [xLimLeft, xLimRight];
        elseif strcmpi(mode,'largestOne')
            diffOut = 0;
            for k = 1:size(all_axes,1)
                for j = 1:size(all_axes,2)
                    if diffOut < diff(all_axes{k,j}.XLim)
                        xLimOut = all_axes{k,j}.XLim;
                        diffOut = diff(all_axes{k,j}.XLim);
                    end
                end
            end
        elseif strcmpi(mode,'smallestOne')
            diffOut = Inf;
            for k = 1:size(all_axes,1)
                for j = 1:size(all_axes,2)
                    if diffOut > diff(all_axes{k,j}.XLim)
                        xLimOut = all_axes{k,j}.XLim;
                        diffOut = diff(all_axes{k,j}.XLim);
                    end
                end
            end            
        else
            error('Wrong mode!');
        end
    end
    
    for k = 1:size(all_axes,1)
        for j = 1:size(all_axes,2)
            set(all_axes{k,j},'XLim',xLimOut);
        end
    end

end