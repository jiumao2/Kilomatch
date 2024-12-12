function [handle_patch, handle_line] = plotShaded(ax, x, y, varargin)
%PLOTSHADED Plot a line with shaded area below
%
% EasyPlot.plotShaded(ax, [2,3,5], [1,0,2;3,4,5]);
% EasyPlot.plotShaded(ax, [1,0,2;3,4,5], [2,3,5], 'axis', 'Y');
% x: 1xn vector if axis == 'X';
% y: 2xn matrix composed of two vectors which defines the area of the patch

alpha = 0.25;
axis = 'X';
lineWidth = EasyPlot.DefaultValue.LineWidth;
lineColor = EasyPlot.DefaultValue.Colororder(1,:);
shadedColor = EasyPlot.DefaultValue.Colororder(1,:);
marker = 'none';
drawLine = 'off';
if nargin>3
    for k = 1:2:size(varargin,2)
        switch varargin{k}
            case 'alpha'
                alpha = varargin{k+1};
            case 'axis'
                axis = varargin{k+1};
            case 'shadedColor'
                shadedColor = varargin{k+1};
            case 'lineColor'
                lineColor = varargin{k+1};
            case 'lineWidth'
                lineWidth = varargin{k+1};
            case 'marker'
                marker = varargin{k+1};
            case 'drawLine'
                drawLine = varargin{k+1};
            otherwise
                error('Unknown argument!')
        end
    end
end

if size(x,1)>size(x,2)
    x=x';
end

if size(y,1)>size(y,2)
    y=y';
end

if size(y,1)==1 && axis == 'X'
    y = [y;zeros(1,length(y))];
elseif size(x,1)==1 && axis == 'Y'
    x = [x;zeros(1,length(x))];
end

if strcmp(drawLine, 'on')
    handle_line = plot(ax, x(1,:),y(1,:),...
        'Color',lineColor,...
        'LineWidth',lineWidth,...
        'Marker',marker);
end

% make closed patch 
if size(y,1)==2
    px = [x, fliplr(x)];
    py = [y(1,:), fliplr(y(2,:))];
elseif size(x,1)==2
    px = [x(1,:), fliplr(x(2,:))];
    py = [y, fliplr(y)]; 
end

% plot shaded area
handle_patch = patch(ax, px, py, 1,...
    'FaceColor', shadedColor,...
    'EdgeColor','none',...
    'FaceAlpha', alpha); 

end