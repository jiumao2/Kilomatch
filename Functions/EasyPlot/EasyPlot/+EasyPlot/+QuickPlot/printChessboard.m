function printChessboard(gridSize)
%PRINTCHESSBOARD Print chessboard picture to A4 paper. The default grid size is 20 mm.

grid_size = 20; % mm
if nargin>=1
    grid_size = gridSize;
end

A4_height = 210;
A4_width = 297;

% rescale to not overflow the screen
scale = 0.1;
grid_size = grid_size*scale;
A4_width = A4_width*scale;
A4_height = A4_height*scale;

num_width = floor(A4_width/grid_size);
num_height = floor(A4_height/grid_size);

if mod(num_width,2) == mod(num_height, 2)
    num_width = num_width-1;
end

margin_horizontal = (A4_width-num_width*grid_size)/2;
margin_vertical = (A4_height-num_height*grid_size)/2;

fig = EasyPlot.figure(...
    'Width', A4_width, ...
    'Height', A4_height, ...
    'MarginBottom', 0, ...
    'MarginLeft', 0, ...
    'MarginRight', 0, ...
    'MarginTop', 0);

fig.Position(1:2) = 0;

ax = EasyPlot.createAxesAgainstFigure(fig, 'leftBottom',...
    'Width', num_width*grid_size,...
    'Height', num_height*grid_size,...
    'MarginLeft', margin_horizontal,...
    'MarginRight', margin_horizontal,...
    'MarginBottom', margin_vertical,...
    'MarginTop', margin_vertical,...
    'XAxisVisible', 'off', ...
    'YAxisVisible', 'off');

img = uint8(zeros(num_height, num_width));
for k = 1:num_height
    for j = 1:num_width
        img(k,j) = uint8(mod(k+j, 2)*255);
    end
end

image(ax, img, ...
    'XData', (1:num_width)-0.5, ...
    'YData', (1:num_height)-0.5);

colormap('gray')
xlim(ax, [0, num_width]);
ylim(ax, [0, num_height]);

set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperOrientation', 'landscape');

EasyPlot.exportFigure(fig, sprintf('Chessboard_%dmm',grid_size/scale), 'dpi', 1200, 'type', 'pdf');
EasyPlot.exportFigure(fig, sprintf('Chessboard_%dmm',grid_size/scale), 'dpi', 1200, 'type', 'png');
close all;
end
