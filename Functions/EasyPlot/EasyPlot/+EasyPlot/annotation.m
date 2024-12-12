function handle = annotation(fig, type, varargin)
    fontName = EasyPlot.DefaultValue.FontName;
    fontSize = EasyPlot.DefaultValue.FontSize;
    fontWeight = EasyPlot.DefaultValue.FontWeight;
    marginLeft = EasyPlot.DefaultValue.AxesMarginLeft;
    marginRight = EasyPlot.DefaultValue.AxesMarginRight;
    marginTop = EasyPlot.DefaultValue.AxesMarginTop;
    marginBottom = EasyPlot.DefaultValue.AxesMarginBottom;

    handle = annotation(fig, type,...
        'Units', EasyPlot.DefaultValue.Units);

    if strcmpi(type, 'textarrow') || strcmpi(type, 'textbox')
        set(handle,...
            'FontName',fontName,...
            'FontWeight',fontWeight,...
            'FontSize', fontSize,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left');
    end

    if strcmpi(type, 'textbox')
        set(handle,...
            'EdgeColor','none',...
            'FitBoxToText','off');        
    end

    handle.UserData.MarginLeft = marginLeft;
    handle.UserData.MarginRight = marginRight;
    handle.UserData.MarginTop = marginTop;
    handle.UserData.MarginBottom = marginBottom;

    if nargin > 2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'MarginLeft')
                handle.UserData.MarginLeft = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginRight')
                handle.UserData.MarginRight = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginTop')
                handle.UserData.MarginTop = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginBottom')
                handle.UserData.MarginBottom = varargin{k+1};
            else
                set(handle,varargin{k},varargin{k+1});
            end
        end
    end

    EasyPlot.Utils.addHandle(fig, handle);
end