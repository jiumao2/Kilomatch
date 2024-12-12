function zLimOut = setZLim(all_axes,mode)
    if nargin<2
        mode = 'largest';
    end
    zLimOut = nan(1,2);

    if ~iscell(all_axes)
        if ~isnumeric(mode) || length(mode)~=2
            error('Wrong mode!')
        else
            zLimOut = mode;
        end
        set(all_axes,'ZLim',zLimOut);
        return
    end

    if isnumeric(mode)
        if length(mode) == 1
            zLimOut = all_axes{mode}.ZLim;
        elseif length(mode) == 2
            zLimOut = mode;
        else
            error('Wrong mode!');
        end
    elseif ischar(mode) || isstring(mode)
        if strcmpi(mode,'largest')
            zLimLeft = Inf;
            zLimRight = -Inf;
            for k = 1:size(all_axes,1)
                for j = 1:size(all_axes,2)
                    left = all_axes{k,j}.ZLim(1);
                    right = all_axes{k,j}.ZLim(2);
                    zLimLeft = min(zLimLeft,left);
                    zLimRight = max(zLimRight,right);
                end
            end
            zLimOut = [zLimLeft, zLimRight];
        elseif strcmpi(mode,'largestOne')
            diffOut = 0;
            for k = 1:size(all_axes,1)
                for j = 1:size(all_axes,2)
                    if diffOut < diff(all_axes{k,j}.ZLim)
                        zLimOut = all_axes{k,j}.ZLim;
                        diffOut = diff(all_axes{k,j}.ZLim);
                    end
                end
            end
        elseif strcmpi(mode,'smallestOne')
            diffOut = Inf;
            for k = 1:size(all_axes,1)
                for j = 1:size(all_axes,2)
                    if diffOut > diff(all_axes{k,j}.ZLim)
                        zLimOut = all_axes{k,j}.ZLim;
                        diffOut = diff(all_axes{k,j}.ZLim);
                    end
                end
            end            
        else
            error('Wrong mode!');
        end
    end
    
    for k = 1:size(all_axes,1)
        for j = 1:size(all_axes,2)
            set(all_axes{k,j},'ZLim',zLimOut);
        end
    end

end