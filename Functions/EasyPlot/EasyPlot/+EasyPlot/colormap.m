function colormap(axes_all, cmap, varargin)
    zeroCenter = 'off';
    zeroPosition = 0;
    if nargin > 2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'zeroCenter')
                zeroCenter = varargin{k+1};
            elseif strcmpi(varargin{k},'zeroPosition')
                zeroPosition = varargin{k+1};
            else
                error('Unknown argument!');
            end
        end
    end
    if ~iscell(axes_all)
        if strcmpi(zeroCenter, 'on')
            cLim = axes_all.CLim - zeroPosition;
            cmap(1:round((2*EasyPlot.DefaultValue.ColormapDivergingWhitePosition-1)*size(cmap,1)),:) = [];
            if cLim(2) > -cLim(1)
                percentCmap = (cLim(2)-cLim(1))./(2*cLim(2));
                cmap = cmap(max(round(size(cmap,1)*(1-percentCmap)), 1):end,:);
            elseif cLim(2) < -cLim(1)
                percentCmap = (cLim(2)-cLim(1))./(2*-cLim(1));
                cmap = cmap(1:min(round(size(cmap,1)*percentCmap), size(cmap, 1)),:);
            end
        end
        colormap(axes_all, cmap);
    else
        % Check whether the CLim are same for all axes
        flag_different = false;
        for k = 1:size(axes_all,1)
            for j = 1:size(axes_all,2)
                if ~all(axes_all{k,j}.CLim == axes_all{1,1}.CLim)
                    warning('The CLim from the input axes are different. It may result to different cmaps for the axes.');
                    flag_different = true;
                    break
                end
            end
            if flag_different
                break
            end
        end
        
        cmap_raw = cmap;
        for k = 1:size(axes_all,1)
            for j = 1:size(axes_all,2)
                if strcmpi(zeroCenter, 'on')
                    cmap = cmap_raw;
                    cLim = axes_all{k,j}.CLim - zeroPosition;
                    cmap(1:round((2*EasyPlot.DefaultValue.ColormapDivergingWhitePosition-1)*size(cmap,1)),:) = [];
                    if cLim(2) > -cLim(1)
                        percentCmap = (cLim(2)-cLim(1))./(2*cLim(2));
                        cmap = cmap(max(round(size(cmap,1)*(1-percentCmap)), 1):end,:);
                    elseif cLim(2) < -cLim(1)
                        percentCmap = (cLim(2)-cLim(1))./(2*-cLim(1));
                        cmap = cmap(1:min(round(size(cmap,1)*percentCmap), size(cmap, 1)),:);
                    end
                end
                colormap(axes_all{k,j}, cmap);
            end
        end
    end  
end
