function handle = setGeneralTitle(axes_all, label, varargin)
    if ~iscell(axes_all)
        axes_all = {axes_all};
    end
    height = 0.8;

    posLeft = 1e8;
    posRight = -1e8;

    fig = axes_all{1}.Parent;
    for k = 1:length(axes_all)
        posLeft = min(posLeft, axes_all{k}.Position(1));
        posRight = max(posRight, axes_all{k}.Position(1)+axes_all{k}.Position(3));
    end
    posTop = axes_all{1}.Position(2)+axes_all{1}.Position(4);

    pos = [posLeft-1,posTop,posRight-posLeft+2,height];
    
    handle = EasyPlot.annotation(fig, 'textbox',...
        'FitBoxToText','on',...
        'String',label,...
        'Position',pos,...
        'EdgeColor','none',...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center',...
        'fontSize',axes_all{1}.XLabel.FontSize,...
        'MarginBottom',0,...
        'MarginLeft',0,...
        'MarginRight',0,...
        'MarginTop',0);

    if nargin>2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'MarginLeft')
                handle.UserData.MarginLeft = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginRight')
                handle.UserData.MarginRight = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginTop')
                handle.UserData.MarginTop = varargin{k+1};
            elseif strcmpi(varargin{k},'MarginBottom')
                handle.UserData.MarginBottom = varargin{k+1};
            elseif strcmpi(varargin{k},'Height')
                handle.Position(4) = varargin{k+1};
            elseif strcmpi(varargin{k},'xShift')
                handle.Position(1) = handle.Position(1)+varargin{k+1};
            elseif strcmpi(varargin{k},'yShift')
                handle.Position(2) = handle.Position(2)+varargin{k+1};
            else
                EasyPlot.set(handle,varargin{k},varargin{k+1});
            end
        end
    end
end