function [ax, ax2, ax3, ax4] = truncAxis(ax, varargin)
if ~isa(ax,'matlab.graphics.axis.Axes')
    ax = gca;
end
fig = ax.Parent;
ax.Units = "normalized";
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';

axisPos = ax.Position;
axisXLim = ax.XLim;
axisYLim = ax.YLim;

axisXScale = diff(axisXLim);
axisYScale = diff(axisYLim);

truncRatio = 1/20;
Xtrunc = [];
Ytrunc = [];
Xratio = NaN; % left/right
Yratio = NaN; % bottom/top
if nargin>1
    for i = 1:2:size(varargin,2)
        switch varargin{i}
            case 'X'
                Xtrunc = varargin{i+1};
            case 'Y'
                Ytrunc = varargin{i+1};
            case 'Xratio'
                Xratio = varargin{i+1};
            case 'Yratio'
                Yratio = varargin{i+1};
            case 'truncRatio'
                truncRatio = varargin{i+1};
            otherwise
                error('wrong argument!');
        end
    end
end

if isempty(Xtrunc)
    % copy axes
    ax2 = copyAxes(ax);

    ax2.XTickLabels = [];
    ax2.XColor = 'none';

    ax.YLim = [axisYLim(1),Ytrunc(1)];
    ax2.YLim = [Ytrunc(2),axisYLim(2)];

    if isnan(Yratio)
        Yratio = (Ytrunc(1)-axisYLim(1))/(axisYLim(2)-Ytrunc(2));
    end

    ax.Position(4) = axisPos(4)*(1-truncRatio)*Yratio/(1+Yratio);

    ax2.Position = ax.Position;
    ax2.Position(2) = axisPos(2)+ax.Position(4)+axisPos(4)*truncRatio;
    ax2.Position(4) = axisPos(4)*(1-truncRatio)*1/(1+Yratio);

    linkaxes([ax,ax2],'x');

    if strcmp(ax.Box,'on')
        ax.Box = 'off';ax2.Box = 'off';
        createLine([1,1].*(ax.Position(1)+ax.Position(3)),[ax.Position(2),ax.Position(2)+ax.Position(4)],'-',ax.Color);
        createLine([1,1].*(ax.Position(1)+ax.Position(3)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],'-',ax.Color);
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax2.Position(2)+ax2.Position(4)),'-',ax.Color);
    else
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax.Position(2)+ax.Position(4)),':',ax.XColor);
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax2.Position(2)),':',ax.XColor);
    end
    createSlash([ax.Position(1)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4]);
    createSlash([ax.Position(1)-.2,ax2.Position(2)-.2,.4,.4]);
%     createSlash([ax.Position(1)+ax.Position(3)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4]);
%     createSlash([ax.Position(1)+ax.Position(3)-.2,ax2.Position(2)-.2,.4,.4]);
    ax.Units = "centimeters";
    ax2.Units = "centimeters";
elseif isempty(Ytrunc) 
    ax2 = copyAxes(ax);

    ax2.YTickLabels = [];
    ax2.YColor = 'none';

    ax.XLim = [axisXLim(1),Xtrunc(1)];
    ax2.XLim = [Xtrunc(2),axisXLim(2)];

    if isnan(Xratio)
        Xratio = (Xtrunc(1)-axisXLim(1))/(axisXLim(2)-Xtrunc(2));
    end

    ax.Position(3) = axisPos(3)*(1-truncRatio)*Xratio/(1+Xratio);

    ax2.Position  =  ax.Position;
    ax2.Position(1) = axisPos(1)+ax.Position(3)+axisPos(3)*truncRatio;
    ax2.Position(3) = axisPos(3)*(1-truncRatio)*1/(1+Xratio);

    linkaxes([ax,ax2],'y')

    if strcmp(ax.Box,'on')
        ax.Box = 'off';ax2.Box = 'off';
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax.Position(2)+ax.Position(4)),'-',ax.XColor);
        createLine([ax2.Position(1),ax2.Position(1)+ax2.Position(3)],[1,1].*(ax.Position(2)+ax.Position(4)),'-',ax.XColor);
        createLine([1,1].*(ax2.Position(1)+ax2.Position(3)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],'-',ax.XColor);
    else
        createLine([1,1].*(ax.Position(1)+ax.Position(3)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],':',ax.XColor);
        createLine([1,1].*(ax2.Position(1)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],':',ax.XColor);
    end
    createSlash([ax.Position(1)+ax.Position(3)-.2,ax.Position(2)-.2,.4,.4]);
    createSlash([ax2.Position(1)-.2,ax.Position(2)-.2,.4,.4]);
%     createSlash([ax.Position(1)+ax.Position(3)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4]);
%     createSlash([ax2.Position(1)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4]);
    ax.Units = "centimeters";
    ax2.Units = "centimeters";
elseif ~isempty(Xtrunc) && ~isempty(Ytrunc)
    ax2 = copyAxes(ax);
    ax3 = copyAxes(ax);
    ax4 = copyAxes(ax);

    ax2.XTickLabels = [];
    ax2.XColor = 'none';
    ax3.XTickLabels = [];
    ax3.XColor = 'none';
    ax3.YTickLabels = [];
    ax3.YColor = 'none';
    ax4.YTickLabels = [];
    ax4.YColor = 'none';

    ax.YLim = [axisYLim(1),Ytrunc(1)];
    ax.XLim = [axisXLim(1),Xtrunc(1)];
    ax2.XLim = [axisXLim(1),Xtrunc(1)];
    ax2.YLim = [Ytrunc(2),axisYLim(2)];
    ax3.XLim = [Xtrunc(2),axisXLim(2)];
    ax3.YLim = [Ytrunc(2),axisYLim(2)];
    ax4.XLim = [Xtrunc(2),axisXLim(2)];
    ax4.YLim = [axisYLim(1),Ytrunc(1)];

    if isnan(Xratio)
        Xratio = (Xtrunc(1)-axisXLim(1))/(axisXLim(2)-Xtrunc(2));
    end    
    if isnan(Yratio)
        Yratio = (Ytrunc(1)-axisYLim(1))/(axisYLim(2)-Ytrunc(2));
    end

    ax.Position(3) = axisPos(3)*(1-truncRatio)*Xratio/(1+Xratio);
    ax.Position(4) = axisPos(4)*(1-truncRatio)*Yratio/(1+Yratio);

    ax2.Position = ax.Position;
    ax2.Position(2) = axisPos(2)+ax.Position(4)+axisPos(4)*truncRatio;
    ax2.Position(3) = axisPos(3)*(1-truncRatio)*Xratio/(1+Xratio);
    ax2.Position(4) = axisPos(4)*(1-truncRatio)*1/(1+Yratio);

    ax3.Position = ax.Position;
    ax3.Position(1) = axisPos(1)+ax.Position(3)+axisPos(3)*truncRatio;
    ax3.Position(2) = axisPos(2)+ax.Position(4)+axisPos(4)*truncRatio;
    ax3.Position(3) = axisPos(3)*(1-truncRatio)*1/(1+Xratio);
    ax3.Position(4) = axisPos(4)*(1-truncRatio)*1/(1+Yratio);

    ax4.Position = ax.Position;
    ax4.Position(1) = axisPos(1)+ax.Position(3)+axisPos(3)*truncRatio;
    ax4.Position(3) = axisPos(3)*(1-truncRatio)*1/(1+Xratio);
    ax4.Position(4) = axisPos(4)*(1-truncRatio)*Yratio/(1+Yratio);

    linkaxes([ax3,ax2],'y');
    linkaxes([ax4,ax3],'x');
    linkaxes([ax,ax2],'x');
    linkaxes([ax,ax4],'y');

    if strcmp(ax.Box,'on')
        ax.Box = 'off';ax2.Box = 'off';ax3.Box = 'off';ax4.Box = 'off';
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax2.Position(2)+ax2.Position(4)),'-',ax.XColor);
        createLine([ax3.Position(1),ax3.Position(1)+ax3.Position(3)],[1,1].*(ax2.Position(2)+ax2.Position(4)),'-',ax.XColor);
        createLine([1,1].*(ax4.Position(1)+ax4.Position(3)),[ax3.Position(2),ax3.Position(2)+ax3.Position(4)],'-',ax.XColor);
        createLine([1,1].*(ax4.Position(1)+ax4.Position(3)),[ax4.Position(2),ax4.Position(2)+ax4.Position(4)],'-',ax.XColor);
    else
        createLine([1,1].*(ax.Position(1)+ax.Position(3)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],':',ax.XColor);
        createLine([1,1].*(ax3.Position(1)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],':',ax.XColor);
        createLine([1,1].*(ax.Position(1)+ax.Position(3)),[ax.Position(2),ax.Position(2)+ax.Position(4)],':',ax.XColor);
        createLine([1,1].*(ax3.Position(1)),[ax.Position(2),ax.Position(2)+ax.Position(4)],':',ax.XColor);
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax.Position(2)+ax.Position(4)),':',ax.XColor);
        createLine([ax.Position(1),ax.Position(1)+ax.Position(3)],[1,1].*(ax2.Position(2)),':',ax.XColor);
        createLine([ax4.Position(1),ax4.Position(1)+ax4.Position(3)],[1,1].*(ax.Position(2)+ax.Position(4)),':',ax.XColor);
        createLine([ax4.Position(1),ax4.Position(1)+ax4.Position(3)],[1,1].*(ax2.Position(2)),':',ax.XColor);
    end
    createSlash([ax.Position(1)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4]);
    createSlash([ax.Position(1)-.2,ax2.Position(2)-.2,.4,.4]);
%     createSlash([ax4.Position(1)+ax4.Position(3)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4]);
%     createSlash([ax4.Position(1)+ax4.Position(3)-.2,ax2.Position(2)-.2,.4,.4]);
    createSlash([ax.Position(1)+ax.Position(3)-.2,ax.Position(2)-.2,.4,.4]);
%     createSlash([ax.Position(1)+ax.Position(3)-.2,ax2.Position(2)+ax2.Position(4)-.2,.4,.4]);
    createSlash([ax4.Position(1)-.2,ax.Position(2)-.2,.4,.4]);
%     createSlash([ax4.Position(1)-.2,ax2.Position(2)+ax2.Position(4)-.2,.4,.4]);

    ax.Units = "centimeters";
    ax2.Units = "centimeters";
    ax3.Units = "centimeters";
    ax4.Units = "centimeters";
end

    function newAX = copyAxes(ax)
        axStruct = get(ax);
        fNames = fieldnames(axStruct);
        newAX = EasyPlot.axes(fig);

        coeList = {'CurrentPoint','XAxis','YAxis','ZAxis','BeingDeleted',...
            'TightInset','NextSeriesIndex','Children','Type','Legend'};
        for n = 1:length(coeList)
            coePos = strcmp(fNames,coeList{n});
            fNames(coePos) = [];
        end
        
        for n = 1:length(fNames)
            newAX.(fNames{n}) = ax.(fNames{n});
        end

        copyobj(ax.Children,newAX)
    end
    
    function handle = createSlash(pos)
        handle = EasyPlot.annotation(fig,'textbox',...
            'String','/',...
            'LineStyle','none',...
            'fontSize',15,...
            'Units','normalized',...
            'Position',pos,...
            'FitBoxToText','on',...
            'VerticalAlignment','middle',...
            'HorizontalAlignment','center',...
            'MarginBottom',0,...
            'MarginRight',0,...
            'MarginLeft',0,...
            'MarginTop',0);

        handle.Units = 'centimeters';
    end

    function handle = createLine(x, y, lineStyle, color)
        handle = EasyPlot.annotation(fig,'line',...
            'LineStyle',lineStyle,...
            'LineWidth',ax.LineWidth,...
            'Color', color,...
            'Units','normalized',...
            'X',x,...
            'Y',y,...
            'MarginBottom',0,...
            'MarginRight',0,...
            'MarginLeft',0,...
            'MarginTop',0);

        handle.Units = 'centimeters';
    end
end