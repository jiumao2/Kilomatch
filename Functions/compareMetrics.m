function compareMetrics(func, matched_pairs, unmatched_pairs, func_names, save_name, varargin)
binwidth = 0.2;
if nargin>5
    for k = 1:2:size(varargin, 2)
        if strcmpi(varargin{k}, 'binwidth')
            binwidth = varargin{k+1};
        end
    end
end

assert(all(size(func) == size(func_names)));

similarity_matched = cell(size(func));
similarity_unmatched = cell(size(func));

for k = 1:size(func, 1)
    for j = 1:size(func, 2)
        similarity_matched{k,j} = zeros(1, size(matched_pairs, 1));
        similarity_unmatched{k,j} = zeros(1, size(unmatched_pairs, 1));
    end
end

disp('Computing the similarity between matched pairs...');
for k = 1:size(matched_pairs, 1)
    for j = 1:size(func, 1)
        for i = 1:size(func, 2)
            similarity_matched{j,i}(k) = func{j,i}(matched_pairs(k,:));
        end
    end
    if mod(k, 1000) == 1
        fprintf('%d / %d done!\n', k, size(matched_pairs, 1));
    end
end

disp('Computing the similarity between unmatched pairs...');
for k = 1:size(unmatched_pairs, 1)
    for j = 1:size(func, 1)
        for i = 1:size(func, 2)
            similarity_unmatched{j,i}(k) = func{j,i}(unmatched_pairs(k,:));
        end
    end
    if mod(k, 1000) == 1
        fprintf('%d / %d done!\n', k, size(unmatched_pairs, 1));
    end
end

% compute the AUC
AUCs = zeros(size(func));
for k = 1:size(func, 1)
    for j = 1:size(func, 2)
        [~, ~, ~, AUCs(k,j)] = perfcurve(...
            [ones(length(similarity_matched{k,j}),1); 2*ones(length(similarity_unmatched{k,j}),1)],...
            [similarity_matched{k,j}'; similarity_unmatched{k,j}'],...
            1);
    end
end

% plot the results
fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, size(func_names, 1), size(func_names, 2),...
    'Width', 3,...
    'Height', 3,...
    'MarginLeft', 1,...
    'MarginTop', 1,...
    'MarginBottom', 1);

for k = 1:size(func, 1)
    for j = 1:size(func, 2)
        histogram(ax_all{k,j}, similarity_unmatched{k,j}, 'FaceColor', 'k', 'BinWidth', binwidth, 'Normalization', 'probability');
        histogram(ax_all{k,j}, similarity_matched{k,j}, 'FaceColor', 'b', 'BinWidth', binwidth, 'Normalization', 'probability');
        title(ax_all{k,j}, [func_names{k,j}, ' (', num2str(AUCs(k,j), '%.3f'), ')']);
    end
end
EasyPlot.setYLim(ax_all);
EasyPlot.setXLim(ax_all);

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, ['./Figures/', save_name]);

end