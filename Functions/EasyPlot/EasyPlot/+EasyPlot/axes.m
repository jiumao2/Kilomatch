function handle = axes(fig, varargin)
    handle = axes(fig,...
        'Units',EasyPlot.DefaultValue.Units,...
        'NextPlot',EasyPlot.DefaultValue.NextPlot,...
        'FontName',EasyPlot.DefaultValue.FontName,...
        'FontSize',EasyPlot.DefaultValue.FontSize);
    handle.UserData.MarginLeft = EasyPlot.DefaultValue.AxesMarginLeft;
    handle.UserData.MarginRight = EasyPlot.DefaultValue.AxesMarginRight;
    handle.UserData.MarginTop = EasyPlot.DefaultValue.AxesMarginTop;
    handle.UserData.MarginBottom = EasyPlot.DefaultValue.AxesMarginBottom;
    handle.Position(3) = EasyPlot.DefaultValue.AxesWidth;
    handle.Position(4) = EasyPlot.DefaultValue.AxesHeight;
    if nargin > 1
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'MarginLeft')
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
            elseif strcmpi(varargin{k},'XAxisVisible')
                handle.XAxis.Visible = varargin{k+1};
            elseif strcmpi(varargin{k},'YAxisVisible')
                handle.YAxis.Visible = varargin{k+1};
            else
                EasyPlot.set(handle,varargin{k},varargin{k+1});
            end
        end
    end
end