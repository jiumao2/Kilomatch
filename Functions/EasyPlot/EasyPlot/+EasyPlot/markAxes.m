function handle = markAxes(fig, axes_all, marker, varargin)
    xShift = 0;
    yShift = 0;
    fontName = EasyPlot.DefaultValue.FontName;
    fontSize = EasyPlot.DefaultValue.FontSizeTitle;
    fontWeight = 'bold';
    marginLeft = EasyPlot.DefaultValue.AxesMarginLeft;
    marginRight = EasyPlot.DefaultValue.AxesMarginRight;
    marginTop = EasyPlot.DefaultValue.AxesMarginTop;
    marginBottom = EasyPlot.DefaultValue.AxesMarginBottom;
    width = 1;
    height = 1;
    horizontalAlignment = 'left';
    verticalAlignment = 'bottom';

    if nargin>3
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'xShift')
                xShift = varargin{k+1};
            elseif strcmpi(varargin{k},'yShift')
                yShift = varargin{k+1};
            elseif strcmpi(varargin{k},'fontSize')
                fontSize = varargin{k+1};
            elseif strcmpi(varargin{k},'fontWeight')
                fontWeight = varargin{k+1};
            elseif strcmpi(varargin{k},'fontName')
                fontName = varargin{k+1};
            elseif strcmpi(varargin{k},'marginLeft')
                marginLeft = varargin{k+1};
            elseif strcmpi(varargin{k},'marginRight')
                marginRight = varargin{k+1};
            elseif strcmpi(varargin{k},'marginTop')
                marginTop = varargin{k+1};
            elseif strcmpi(varargin{k},'marginBottom')
                marginBottom = varargin{k+1};
            elseif strcmpi(varargin{k},'Width')
                width = varargin{k+1};
            elseif strcmpi(varargin{k},'Height')
                height = varargin{k+1};
            elseif strcmpi(varargin{k},'HorizontalAlignment')
                horizontalAlignment = varargin{k+1};
            elseif strcmpi(varargin{k},'VerticalAlignment')
                verticalAlignment = varargin{k+1};
            else
                error('Wrong argument!');
            end
        end
    end

    if ~iscell(axes_all)
        pos = zeros(1,4);
        pos(3) = width;
        pos(4) = height;
        pos(1) = axes_all.Position(1) - pos(3) + xShift;
        pos(2) = axes_all.Position(2) + axes_all.Position(4) + yShift;
        handle = EasyPlot.annotation(fig,'textbox',...
            'Position',pos,...
            'String',marker,...
            'FontName',fontName,...
            'FontWeight',fontWeight,...
            'FontSize',fontSize,...
            'HorizontalAlignment', horizontalAlignment,...
            'VerticalAlignment', verticalAlignment);

        handle.UserData.MarginLeft = marginLeft;
        handle.UserData.MarginRight = marginRight;
        handle.UserData.MarginTop = marginTop;
        handle.UserData.MarginBottom = marginBottom;
        return
    end

    for k = 1:length(axes_all)
        pos = zeros(1,4);
        pos(3) = width;
        pos(4) = height;
        pos(1) = axes_all{k}.Position(1) - pos(3) + xShift;
        pos(2) = axes_all{k}.Position(2) + axes_all{k}.Position(4) + yShift;   
        handle(k) = EasyPlot.annotation(fig,'textbox',...
            'Position',pos,...
            'String',marker{k},...
            'FontName',fontName,...
            'FontWeight',fontWeight,...
            'FontSize', fontSize,...
            'HorizontalAlignment', horizontalAlignment,...
            'VerticalAlignment', verticalAlignment);

        handle(k).UserData.MarginLeft = marginLeft;
        handle(k).UserData.MarginRight = marginRight;
        handle(k).UserData.MarginTop = marginTop;
        handle(k).UserData.MarginBottom = marginBottom;
    end
end