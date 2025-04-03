function handle_all = scalebar(ax, type, varargin)
%SCALEBAR creates a scalebar on an axes
%
% SCALEBAR(ax)
% SCALEBAR(ax, type)
% SCALEBAR(ax, type, name, value)
% handle_all = SCALEBAR(___)
% 
% Draws X / Y scalebars on the axes and returns the handles to the scalebars.
% The Unit property of the axes must be set to 'centimeters'.
% 
% INPUT
%     ax:       handle to the axes
%     type:     the type of scalebar
%               'X' (default) | 'Y' | 'XY'
%
% NAME-VALUE pairs
% 
%     Location: location of the scalebar. Possible values are
%               southwest (default)
%               southeast
%               northeast
%               northwest
% 
%     xBarLabel: the label of the X scalebar (e.g., '1 mm')
%     yBarLabel: the label of the Y scalebar
% 
%     xBarLength: the length of the X scalebar (e.g., for xBarLabel '2 mm', the xBarLength is 2)
%     yBarLength: the length of the Y scalebar
% 
%     xBarRatio: the ratio between unit of xBarLength and the unit of the X axis (default, 1)
%         e.g., the unit of X axis is cm and the unit of X scalebar is mm, then the xBarRatio should be 0.1
%         
%     yBarRatio: the ratio between unit of yBarLength and the unit of the Y axis (default, 1)
% 
%     xTextDistance: to increase or decrease the distance between text and X scalebar (default, 0)
%     yTextDistance: to increase or decrease the distance between text and Y scalebar (default, 0)
% 
%     xTextPosition: location of the X text, on the top / bottom of the X scalebar (default, best fit to the scalebar) 
%             top | bottom
%     yTextPosition: location of the Y text, on the left / right of the Y scalebar (default, best fit to the scalebar) 
%             left | right
% 
%     showXText: whether to show the xBarLabel
%             'on' (default) | 'off'
%     showYText: whether to show the xBarLabel
%             'on' (default) | 'off'
% 
% 
%     color: color of scalebar (default, black)
% 
%     lineWidth: the lineWidth of the scalebar (default, 2)
% 
%     fontWeight: the fontWeight of the text
%             'bold' | 'normal' (default)
% 
%     fontSize: the fontSize of the text (default, 8)
% 
%     fontName = the fontName of the text (default, 'Arial')
% 
% 
% OUTPUT
%     handle_all: all handles genenrated by the funcion in 1xn cell
%
%
% Created 14 Jan 2024 by Yue Huang
%

color = 'k';
lineWidth = 2;
fontWeight = 'normal';
fontSize = EasyPlot.DefaultValue.FontSize;
fontName = EasyPlot.DefaultValue.FontName;
xBarLabel = '';
yBarLabel = '';
xBarRatio = 1; % realUnit / xUnit
yBarRatio = 1; % realUnit / yUnit
xBarLength = 1;
yBarLength = 1;
showXText = 'on';
showYText = 'on';
location = 'southwest';

xBarDirection = '';
yBarDirection = '';
xTextPosition = ''; % top or bottom
yTextPosition = ''; % left or right
xTextDistance = 0;
yTextDistance = 0;

if nargin < 2
    type = 'X';
end

if nargin > 2
    for k = 1:2:size(varargin,2)
        if strcmpi(varargin{k},'color')
            color = varargin{k+1};
        elseif strcmpi(varargin{k},'lineWidth')
            lineWidth = varargin{k+1};
        elseif strcmpi(varargin{k},'fontWeight')
            fontWeight = varargin{k+1};
        elseif strcmpi(varargin{k},'fontSize')
            fontSize = varargin{k+1};
        elseif strcmpi(varargin{k},'xBarLabel')
            xBarLabel = varargin{k+1};
        elseif strcmpi(varargin{k},'yBarLabel')
            yBarLabel = varargin{k+1};
        elseif strcmpi(varargin{k},'xBarRatio')
            xBarRatio = varargin{k+1};
        elseif strcmpi(varargin{k},'yBarRatio')
            yBarRatio = varargin{k+1};
        elseif strcmpi(varargin{k},'xBarLength')
            xBarLength = varargin{k+1};
        elseif strcmpi(varargin{k},'yBarLength')
            yBarLength = varargin{k+1};
        elseif strcmpi(varargin{k},'showXText')
            showXText = varargin{k+1};
        elseif strcmpi(varargin{k},'showYText')
            showYText = varargin{k+1};
        elseif strcmpi(varargin{k},'location')
            location = varargin{k+1};
        elseif strcmpi(varargin{k},'xTextPosition')
            xTextPosition = varargin{k+1};
        elseif strcmpi(varargin{k},'yTextPosition')
            yTextPosition = varargin{k+1};
        elseif strcmpi(varargin{k},'xTextDistance')
            xTextDistance = varargin{k+1};
        elseif strcmpi(varargin{k},'yTextDistance')
            yTextDistance = varargin{k+1};
        else
            error('Wrong argument!');
        end
    end
end

% deal with textPositio
if isempty(xTextPosition)
    if contains(location, 'north')
        xTextPosition = 'top';
    else
        xTextPosition = 'bottom';
    end
end

if isempty(yTextPosition)
    if contains(location, 'west')
        yTextPosition = 'left';
    else
        yTextPosition = 'right';
    end
end

% deal with barDirection
if isempty(xBarDirection)
    if contains(location, 'west')
        xBarDirection = 'right';
    else
        xBarDirection = 'left';
    end
end 

if isempty(yBarDirection)
    if contains(location, 'north')
        yBarDirection = 'bottom';
    else
        yBarDirection = 'top';
    end
end 

% get the anchor
x_pos_range = [ax.Position(1), ax.Position(1)+ax.Position(3)];
y_pos_range = [ax.Position(2), ax.Position(2)+ax.Position(4)];
switch location
    case 'northeast'
        anchor = [x_pos_range(2) - range(x_pos_range)*0.05, ...
                  y_pos_range(2) - range(y_pos_range)*0.05];
    case 'northwest'
        anchor = [x_pos_range(1) + range(x_pos_range)*0.05, ...
                  y_pos_range(2) - range(y_pos_range)*0.05];
    case 'southwest'
        anchor = [x_pos_range(1) + range(x_pos_range)*0.05, ...
                  y_pos_range(1) + range(y_pos_range)*0.05];
    case 'southeast'
        anchor = [x_pos_range(2) - range(x_pos_range)*0.05, ...
                  y_pos_range(1) + range(y_pos_range)*0.05];
end   

% calculate the length of line
xLineLength = ax.Position(3)./range(ax.XLim)*xBarLength*xBarRatio;
yLineLength = ax.Position(4)./range(ax.YLim)*yBarLength*yBarRatio;

% calculate the positions of line and text
if strcmpi(xBarDirection, 'right')
    x_xLine = [anchor(1), anchor(1)+xLineLength];
    y_xLine = [anchor(2), anchor(2)];
else
    x_xLine = [anchor(1)-xLineLength, anchor(1)];
    y_xLine = [anchor(2), anchor(2)];
end

if strcmpi(yBarDirection, 'top')
    x_yLine = [anchor(1), anchor(1)];
    y_yLine = [anchor(2), anchor(2)+yLineLength];
else
    x_yLine = [anchor(1), anchor(1)];
    y_yLine = [anchor(2)-yLineLength, anchor(2)];
end

height_textbox = 1.5;
width_textbox = 2;
if strcmpi(xTextPosition, 'bottom')
    xTextboxPosition = [...
        x_xLine(1)+0.5*xLineLength-0.5*width_textbox,...
        y_xLine(1)-height_textbox-xTextDistance,...
        width_textbox,...
        height_textbox];
else
    xTextboxPosition = [...
        x_xLine(1)+0.5*xLineLength-0.5*width_textbox,...
        y_xLine(1)+xTextDistance,...
        width_textbox,...
        height_textbox];
end

if strcmpi(yTextPosition, 'left')
    yTextboxPosition = [...
        x_yLine(1)-width_textbox-yTextDistance,...
        y_yLine(1)+0.5*yLineLength-0.5*height_textbox,...
        width_textbox,...
        height_textbox];
else
    yTextboxPosition = [...
        x_yLine(1)+yTextDistance,...
        y_yLine(1)+0.5*yLineLength-0.5*height_textbox,...
        width_textbox,...
        height_textbox];
end


% draw scalebar
fig = ax.Parent;

handle_all = {};
if contains(type, 'X')
    h_xLine = EasyPlot.annotation(fig, 'line',...
        'X', x_xLine,...
        'Y', y_xLine,...
        'LineWidth', lineWidth,...
        'Color', color,...
        'MarginLeft', 0,...
        'MarginRight', 0,...
        'MarginTop', 0,...
        'MarginBottom', 0);
    handle_all = [handle_all, {h_xLine}];

    if strcmpi(showXText, 'on')
        if strcmpi(xTextPosition, 'bottom')
            x_verticalAlignment = 'top';
        else
            x_verticalAlignment = 'bottom';
        end
        h_xText = EasyPlot.annotation(fig, 'textbox',...
            'String', xBarLabel,...
            'fontName', fontName,...
            'fontSize', fontSize,...
            'fontWeight', fontWeight,...
            'Color', color,...
            'Position', xTextboxPosition,...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', x_verticalAlignment,...
            'MarginLeft', 0,...
            'MarginRight', 0,...
            'MarginTop', 0,...
            'MarginBottom', 0,...
            'FitBoxToText', 'on'); % make the width and height fit to the textbox
        handle_all = [handle_all, {h_xText}];
    end
end

if contains(type, 'Y')
    h_yLine = EasyPlot.annotation(fig, 'line',...
        'X', x_yLine,...
        'Y', y_yLine,...
        'LineWidth', lineWidth,...
        'Color', color,...
        'MarginLeft', 0,...
        'MarginRight', 0,...
        'MarginTop', 0,...
        'MarginBottom', 0);
    handle_all = [handle_all, {h_yLine}];

    if strcmpi(showYText, 'on')
        if strcmpi(yTextPosition, 'left')
            y_horizontalAlignment = 'right';
        else
            y_horizontalAlignment = 'left';
        end
        h_yText = EasyPlot.annotation(fig, 'textbox',...
            'String', yBarLabel,...
            'fontName', fontName,...
            'fontSize', fontSize,...
            'fontWeight', fontWeight,...
            'Color', color,...
            'Position', yTextboxPosition,...
            'HorizontalAlignment', y_horizontalAlignment,...
            'VerticalAlignment', 'middle',...
            'MarginLeft', 0,...
            'MarginRight', 0,...
            'MarginTop', 0,...
            'MarginBottom', 0,...
            'FitBoxToText', 'on'); % make the width and height fit to the textbox
        handle_all = [handle_all, {h_yText}];
    end
end

end