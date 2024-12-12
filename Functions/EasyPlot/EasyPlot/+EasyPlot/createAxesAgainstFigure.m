function handle = createAxesAgainstFigure(fig, location, varargin)
    handle = EasyPlot.axes(fig);
    if nargin > 3
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

    if strcmpi(location,'leftTop')
        handle.Position(1) = fig.UserData.MarginLeft+handle.UserData.MarginLeft;
        handle.Position(2) = fig.Position(4)-fig.UserData.MarginTop-handle.UserData.MarginTop-handle.Position(4);
    elseif strcmpi(location,'rightTop')
        handle.Position(1) = fig.Position(3)-fig.UserData.MarginRight-handle.UserData.MarginRight-handle.Position(3);
        handle.Position(2) = fig.Position(4)-fig.UserData.MarginTop-handle.UserData.MarginTop-handle.Position(4);
    elseif strcmpi(location, 'leftBottom')
        handle.Position(1) = fig.UserData.MarginLeft+handle.UserData.MarginLeft;
        handle.Position(2) = fig.UserData.MarginBottom+handle.UserData.MarginBottom;
    elseif strcmpi(location,'rightBottom')
        handle.Position(1) = fig.Position(3)-fig.UserData.MarginRight-handle.UserData.MarginRight-handle.Position(3);
        handle.Position(2) = fig.UserData.MarginBottom+handle.UserData.MarginBottom;
    end
end