function axes_all_copied = copyAxes(axes_all, location, varargin)
%COPYAXES Copy the position and margins of existing axes. But not all the
%properties are copied.
    if ~iscell(axes_all)
        fig = axes_all.Parent;
        axes_all_copied = EasyPlot.createAxesAgainstAxes(fig, axes_all, location, varargin{:});
        return
    end
    fig = axes_all{1,1}.Parent;

    top = -Inf;
    bottom = Inf;
    left = Inf;
    right = -Inf;
    for k = 1:size(axes_all, 1)
        for j = 1:size(axes_all, 2)
            top_this = axes_all{k,j}.Position(2)+axes_all{k,j}.Position(4)+axes_all{k,j}.UserData.MarginTop;
            bottom_this = axes_all{k,j}.Position(2)-axes_all{k,j}.UserData.MarginBottom;
            left_this = axes_all{k,j}.Position(1)-axes_all{k,j}.UserData.MarginLeft;
            right_this = axes_all{k,j}.Position(1)+axes_all{k,j}.Position(3)+axes_all{k,j}.UserData.MarginRight;
            
            top = max(top,top_this);
            bottom = min(bottom,bottom_this);
            left = min(left,left_this);
            right = max(right,right_this);
        end
    end

    if strcmpi(location, 'left')
        dx = left-right;
        dy = 0;
    elseif strcmpi(location, 'right')
        dx = right-left;
        dy = 0;
    elseif strcmpi(location, 'top')
        dx = 0;
        dy = top-bottom;
    elseif strcmpi(location, 'bottom')
        dx = 0;
        dy = bottom-top;
    elseif strcmpi(location, 'leftTop')
        dx = left-right;
        dy = top-bottom;
    elseif strcmpi(location, 'rightTop')
        dx = right-left;
        dy = top-bottom;
    elseif strcmpi(location, 'leftBottom')
        dx = left-right;
        dy = bottom-top;
    elseif strcmpi(location, 'rightBottom')
        dx = right-left;
        dy = bottom-top;
    else
        error('Wrong location!');
    end
    
    axes_all_copied = cell(size(axes_all));
    for k = 1:size(axes_all, 1)
        for j = 1:size(axes_all, 2)
            ax_this = axes_all{k,j};

            axStruct = get(ax_this);
            fNames = fieldnames(axStruct);
            ax_new = EasyPlot.axes(fig);
    
            coeList = {'CurrentPoint','XAxis','YAxis','ZAxis','BeingDeleted',...
            'TightInset','NextSeriesIndex','Children','Type','Legend',...
            'XLabel','YLabel','ZLabel','Title'...
            'XLim','YLim','ZLim','CLim','Alim'};
            for i = 1:length(coeList)
                coePos = strcmp(fNames,coeList{i});
                fNames(coePos) = [];
            end
            
            for i = 1:length(fNames)
                ax_new.(fNames{i}) = ax_this.(fNames{i});
            end
    
            axes_all_copied{k,j} = ax_new;
        end
    end

    EasyPlot.move(axes_all_copied, 'dx', dx, 'dy', dy);
    EasyPlot.set(axes_all_copied,...
        'XLimMode', 'auto',...
        'YLimMode', 'auto',...
        'ZLimMode', 'auto',...
        'CLimMode', 'auto');
end