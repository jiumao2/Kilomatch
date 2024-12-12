function handle = colorbar(axes, varargin)
    % keep the position of the axes unchanged
    position_raw = axes.Position;

    handle = colorbar(axes,...
        'Units', EasyPlot.DefaultValue.Units,...
        'FontName',EasyPlot.DefaultValue.FontName,...
        'FontSize',EasyPlot.DefaultValue.FontSize);

    position_raw_colorbar = handle.Position;
    handle.Location = 'Manual';
    handle.Position(1) = position_raw_colorbar(1)+position_raw(3)-axes.Position(3);
    handle.Position(2) = position_raw_colorbar(2)+position_raw(4)-axes.Position(4);
    handle.Position(3:4) = position_raw_colorbar(3:4);
    axes.Position = position_raw;

    handle.UserData.MarginLeft = EasyPlot.DefaultValue.AxesMarginLeft;
    handle.UserData.MarginRight = EasyPlot.DefaultValue.AxesMarginRight;
    handle.UserData.MarginTop = EasyPlot.DefaultValue.AxesMarginTop;
    handle.UserData.MarginBottom = EasyPlot.DefaultValue.AxesMarginBottom;

    cmap = axes.Colormap;
    zeroCenter = 'off';
    zeroPosition = 0;

    if nargin>1
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k}, 'label')
                ylabel(handle,varargin{k+1});
            elseif strcmpi(varargin{k}, 'colormap')
                cmap = varargin{k+1};
            elseif strcmpi(varargin{k}, 'zeroCenter')
                zeroCenter = varargin{k+1};
            elseif strcmpi(varargin{k}, 'zeroPosition')
                zeroPosition = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginLeft')
                handle.UserData.MarginLeft = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginRight')
                handle.UserData.MarginRight = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginTop')
                handle.UserData.MarginTop = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginBottom')
                handle.UserData.MarginBottom = varargin{k+1};
            elseif strcmpi(varargin{k},'Width')
                handle.Position(3) = varargin{k+1};
            elseif strcmpi(varargin{k},'Height')
                handle.Position(4) = varargin{k+1};
            elseif strcmpi(varargin{k},'Ticks')
                ticks = sort(varargin{k+1});
                EasyPlot.set(handle, 'Ticks', ticks);
            elseif strcmpi(varargin{k},'TickLabels')
                labels = varargin{k+1};
                if isnumeric(labels)
                    labels_out = cell(1, length(labels));
                    for i_label = 1:length(labels_out)
                        labels_out{i_label} = num2str(labels(i_label));
                    end
                elseif iscell(labels) && isnumeric(labels{1})
                    labels_out = cell(1, length(labels));
                    for i_label = 1:length(labels_out)
                        labels_out{i_label} = num2str(labels{i_label});
                    end
                else
                    labels_out = labels;
                end
                EasyPlot.set(handle, 'TickLabels', labels_out);
            else
                EasyPlot.set(handle,varargin{k},varargin{k+1});
            end
        end
    end
    
    EasyPlot.colormap(axes, cmap, 'zeroCenter', zeroCenter, 'zeroPosition', zeroPosition);
end