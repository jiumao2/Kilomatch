function handle = legend(ax, legendNames, varargin)
    fontName = EasyPlot.DefaultValue.FontName;
    fontSize = EasyPlot.DefaultValue.FontSize;
    fontWeight = EasyPlot.DefaultValue.FontWeight;
    marginLeft = 0;
    marginRight = 0;
    marginTop = 0;
    marginBottom = 0;
    location = 'NorthEast';
    selectedPlots = [];
    Box = 'off';
    lineLength = NaN; % the scale to original lineLength
    lineWidth = NaN;
    markerSize = NaN;

    if nargin > 2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'Location')
                location = varargin{k+1};
            elseif strcmpi(varargin{k},'selectedPlots')
                selectedPlots = varargin{k+1};
            elseif strcmpi(varargin{k},'Box')
                Box = varargin{k+1};
            elseif strcmpi(varargin{k},'lineLength')
                lineLength = varargin{k+1};
            elseif strcmpi(varargin{k},'lineWidth')
                lineWidth = varargin{k+1};
            elseif strcmpi(varargin{k},'markerSize')
                markerSize = varargin{k+1};
            end
        end
    end    
    
    ax_position_raw = ax.Position;
    width_raw = ax.Position(3);
    height_raw = ax.Position(4);
    if isempty(selectedPlots)   
        [handle, labelhandles] = legend(ax, legendNames,...
            'Location', location,...
            'Box', Box,...
            'Units', EasyPlot.DefaultValue.Units,...
            'FontName', fontName,...
            'FontSize', fontSize,...
            'FontWeight', fontWeight);
    else
        [handle, labelhandles] = legend(ax, selectedPlots, legendNames,...
            'Location', location,...
            'Box', Box,...
            'Units', EasyPlot.DefaultValue.Units,...
            'FontName', fontName,...
            'FontSize', fontSize,...
            'FontWeight', fontWeight);
    end

%     handle.ItemTokenSize(1) = lineLength;
    for k = 1:length(labelhandles)
        if strcmpi(labelhandles(k).Type, 'Line')
            if length(labelhandles(k).XData) ~= 2
                continue
            end
            
            if ~isnan(lineLength)
                labelhandles(k).XData(1) = labelhandles(k).XData(1) + (1-lineLength)*range(labelhandles(k).XData);
            end
            if ~isnan(lineWidth)
                labelhandles(k).LineWidth = lineWidth;
            end
            if ~isnan(markerSize)
                labelhandles(k).MarkerSize = markerSize;
            end
        elseif strcmpi(labelhandles(k).Type, 'Text')
            labelhandles(k).FontName = fontName;
            labelhandles(k).FontSize = fontSize;
            labelhandles(k).FontWeight = fontWeight;
        end
    end

    if nargin > 2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'MarginLeft')
                marginLeft = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginRight')
                marginRight = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginTop')
                marginTop = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginBottom')
                marginBottom = varargin{k+1};
            elseif strcmpi(varargin{k},'Location')
                
            elseif strcmpi(varargin{k},'selectedPlots')

            elseif strcmpi(varargin{k},'Box')

            elseif strcmpi(varargin{k},'lineLength')

            else
                EasyPlot.set(handle,varargin{k},varargin{k+1});
            end
        end
    end
    
    position_raw = handle.Position;

    handle.Location = 'None';
    handle.Position = position_raw;
    dx = width_raw - ax.Position(3);
    dy = height_raw - ax.Position(4);
    ax.Position = ax_position_raw;
    EasyPlot.move(handle, 'dx', dx, 'dy', dy);

    handle.UserData.MarginLeft = marginLeft;
    handle.UserData.MarginRight = marginRight;
    handle.UserData.MarginTop = marginTop;
    handle.UserData.MarginBottom = marginBottom;

end