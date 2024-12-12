function handle_all = significanceLine(ax, xStart, xEnd, height, varargin)
%SIGNIFICANCELINE creates a line that marks the statistical significance 
%
% SIGNIFICANCELINE(ax, xStart, xEnd, height)
% SIGNIFICANCELINE(ax, xStart, xEnd, height, name, value)
% handle_all = SCALEBAR(___)
% 
% Draws a line that marks the statistical significance on the axes and 
% returns the handles to the lines and text.
% 
% INPUT
%     ax:       handle to the axes
%     xStart:   the start position of the line
%     xEnd:     the start position of the line
%     height:   the Y position of the line
%
% NAME-VALUE pairs
% 
%     marker: the text shown on the line (default, '*')
% 
%     edgeLength: the length of the line on the edges.
%             The edge will not shown if edgeLength is 0 (default, 0.1*line_length)
% 
%     markerSpace: the space between the marker and the line (default, 0)
% 
%     color: color of the line and text (default, black)
% 
%     lineWidth: the lineWidth of the significance line (default, 1)
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
% Created 23 Jan 2024 by Yue Huang
%

marker = '*';
edgeLength = 0.1*(xEnd - xStart);
markerSpace = 0;
lineWidth = 1;
fontWeight = 'normal';
fontSize = EasyPlot.DefaultValue.FontSize;
fontName = EasyPlot.DefaultValue.FontName;
color = 'black';

if nargin > 2
    for k = 1:2:size(varargin,2)
        if strcmpi(varargin{k},'color')
            color = varargin{k+1};
        elseif strcmpi(varargin{k},'lineWidth')
            lineWidth = varargin{k+1};
        elseif strcmpi(varargin{k},'fontName')
            fontName = varargin{k+1};    
        elseif strcmpi(varargin{k},'fontWeight')
            fontWeight = varargin{k+1};
        elseif strcmpi(varargin{k},'fontSize')
            fontSize = varargin{k+1};
        elseif strcmpi(varargin{k},'marker')
            marker = varargin{k+1};
        elseif strcmpi(varargin{k},'edgeLength')
            edgeLength = varargin{k+1};
        elseif strcmpi(varargin{k},'markerSpace')
            markerSpace = varargin{k+1};
        else
            error('Wrong argument!');
        end
    end
end

% draw the line and text
handle_line = plot(ax, [xStart, xEnd], height*ones(1,2), '-',...
    'Color', color,...
    'LineWidth', lineWidth);
handle_text = text(ax, 0.5*(xStart+xEnd), height+markerSpace, marker,...
    'FontName', fontName,...
    'FontSize', fontSize,...
    'FontWeight', fontWeight,...
    'Color', color,...
    'HorizontalAlignment', 'center');


if edgeLength == 0
    handle_all = {handle_line, handle_text};
else
    handle_left_edge = plot(ax, [xStart, xStart], [height-edgeLength, height], '-',...
        'Color', color,...
        'LineWidth', lineWidth);
    handle_right_edge = plot(ax, [xEnd, xEnd], [height-edgeLength, height], '-',...
        'Color', color,...
        'LineWidth', lineWidth);

    handle_all = {handle_line, handle_left_edge, handle_right_edge, handle_text};
end

end