function axes_all = createGridAxes(fig, nRow, nColumn, varargin)
marginLeft = EasyPlot.DefaultValue.AxesMarginLeft;
marginRight = EasyPlot.DefaultValue.AxesMarginRight;
marginTop = EasyPlot.DefaultValue.AxesMarginTop;
marginBottom = EasyPlot.DefaultValue.AxesMarginBottom;
XAxisVisible = 'on';
YAxisVisible = 'on';

width = 3;
height = 3;
otherNames = {};
otherValues = {};
if nargin > 3
    for k = 1:2:size(varargin,2)
        switch varargin{k}
            case 'Width'
                width = varargin{k+1};
            case 'Height'
                height = varargin{k+1};
            case 'MarginLeft'
                marginLeft = varargin{k+1};
            case 'MarginRight'
                marginRight = varargin{k+1};
            case 'MarginTop'
                marginTop = varargin{k+1};
            case 'MarginBottom'
                marginBottom = varargin{k+1};
            case 'XAxisVisible'
                XAxisVisible = varargin{k+1};
            case 'YAxisVisible'
                YAxisVisible = varargin{k+1};
            otherwise
                otherNames{end+1} = varargin{k};
                otherValues{end+1} = varargin{k+1};
        end
    end
end

axes_all = cell(nRow,nColumn);
for k = 1:nRow
    for j = 1:nColumn
        position = [...
            fig.UserData.MarginLeft+marginLeft+(j-1)*(width+marginLeft+marginRight),...
            fig.UserData.MarginBottom+marginBottom+(nRow-k)*(height+marginTop+marginBottom),...
            width,...
            height];
        axes_all{k,j} = EasyPlot.axes(fig,...
            "Position",position,...
            "MarginTop",marginTop,...
            "MarginBottom",marginBottom,...
            "MarginLeft",marginLeft,...
            "MarginRight",marginRight,...
            'XAxisVisible',XAxisVisible,...
            'YAxisVisible',YAxisVisible);
    end
end

for k = 1:length(otherNames)
    EasyPlot.set(axes_all, otherNames{k}, otherValues{k});
end
end