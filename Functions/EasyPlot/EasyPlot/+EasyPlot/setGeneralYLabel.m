function handle = setGeneralYLabel(axes_all, label, varargin)
    if ~iscell(axes_all)
        axes_all = {axes_all};
    end
    width = 1.5;
    height = 5;
    if nargin>2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'Height')
                height = varargin{k+1};
            elseif strcmpi(varargin{k},'Width')
                width = varargin{k+1};
            end
        end
    end

    posTop = -1e8;
    posBottom = 1e8;

    fig = axes_all{1}.Parent;
    for k = 1:length(axes_all)
        posTop = max(posTop, axes_all{k}.Position(2)+axes_all{k}.Position(4));
        posBottom = min(posBottom, axes_all{k}.Position(2));
    end
    posLeft = axes_all{1}.Position(1)-width;
    pos_vertical_center = 0.5*(posTop+posBottom);

    % compute the position after rotating 90 degree at left-bottom corner
    pos = [posLeft+width, pos_vertical_center-height/2, height, width];
    
    handle = EasyPlot.annotation(fig, 'textbox',...
        'FitBoxToText','on',...
        'String',label,...
        'Position',pos,...
        'EdgeColor','none',...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center',...
        'fontSize',axes_all{1}.YLabel.FontSize,...
        'MarginBottom',0,...
        'MarginLeft',0,...
        'MarginRight',0,...
        'MarginTop',0,...
        'Rotation',90);

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
            elseif strcmpi(varargin{k},'xShift')
                handle.Position(1) = handle.Position(1)+varargin{k+1}; % because of rotation
            elseif strcmpi(varargin{k},'yShift')
                handle.Position(2) = handle.Position(2)+varargin{k+1}; % because of rotation
            elseif strcmpi(varargin{k},'Height')

            elseif strcmpi(varargin{k},'Width')

            else
                EasyPlot.set(handle,varargin{k},varargin{k+1});
            end
        end
    end
end