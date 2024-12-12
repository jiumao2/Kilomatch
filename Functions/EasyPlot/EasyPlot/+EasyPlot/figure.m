function handle = figure(varargin)
    handle = figure(...
        "Units",EasyPlot.DefaultValue.Units,...
        "NextPlot",EasyPlot.DefaultValue.NextPlot,...
        'Renderer',EasyPlot.DefaultValue.Renderer);
    handle.UserData.MarginLeft = EasyPlot.DefaultValue.FigureMarginLeft;
    handle.UserData.MarginRight = EasyPlot.DefaultValue.FigureMarginRight;
    handle.UserData.MarginTop = EasyPlot.DefaultValue.FigureMarginTop;
    handle.UserData.MarginBottom = EasyPlot.DefaultValue.FigureMarginBottom;
    if nargin > 0
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
            else
                set(handle,varargin{k},varargin{k+1});
            end
        end
    end
end