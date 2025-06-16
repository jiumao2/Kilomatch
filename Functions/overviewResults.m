function overviewResults(user_settings, Output)
% (1) Unit number across sessions (position on the probe)
% (2) Motion and matched units' position
% (3) Similairty distribution
% (4) 2D scatter of similarity
% (5) Clustering result (unsorted, sorted)
% (6) P(matched) vs. delta session
% (7) Presence of unique neurons

n_session = Output.NumSession;
n_cluster = Output.NumClusters;
sessions = Output.Sessions;
locations = Output.Locations(:, 1:2);

fig = EasyPlot.figure();

% Unit number across sessions
ax_session = EasyPlot.axes(fig,...
    'Width', 6,...
    'Height', 3,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 1);

x_min = min(locations(:,1));
x_max = max(locations(:,1));
x_scale = 0.8;

x_plot = [];
y_plot = [];
for k = 1:length(sessions)
    x_this = (locations(k,1) - x_min)./(x_max-x_min);
    x_this = (x_this-0.5)*x_scale;
    x_plot = [x_plot, x_this + sessions(k)];
    y_plot = [y_plot, locations(k,2)];
end

n_units_each_session = zeros(1, n_session);
for k = 1:n_session
    n_units_each_session(k) = sum(sessions == k);
end

color_axis = [0.5,0.5,0.5];
yyaxis(ax_session, 'right');
plot(ax_session, x_plot, y_plot, '.', 'Color', color_axis);
ax_session.YAxis(2).Color = color_axis;
ylabel(ax_session, 'Depth (um)');

yyaxis(ax_session, 'left');
plot(ax_session, 1:n_session, n_units_each_session, 'k-', 'LineWidth', 2);
ax_session.YAxis(1).Color = [0,0,0];
ylabel(ax_session, 'Number of units');

xlim(ax_session, [0.5, n_session+0.5]);
xlabel(ax_session, 'Sessions');

% Motion
probe_positions = median(Output.Motion.LinearScale*Output.Motion.Linear*mean(locations(:,2)) + Output.Motion.Constant, 1);

x_plot = [];
y_plot = [];
for k = 1:n_cluster
    units = find(Output.IdxCluster == k);
    sessions_this = sessions(units);

    x_this = (locations(units,1)' - x_min)./(x_max-x_min);
    x_this = (x_this-0.5)*x_scale;

    depth = locations(units, 2)';
    depth_plot = depth - mean(depth) + mean(probe_positions(sessions_this));
    x_plot = [x_plot, x_this + sessions_this];
    y_plot = [y_plot, depth_plot];
end

ax_motion = EasyPlot.createAxesAgainstAxes(fig, ax_session, 'bottom',...
    'MarginBottom', 1);
plot(ax_motion, x_plot, y_plot, '.', 'Color', [0.5,0.5,0.5]);
plot(ax_motion, 1:n_session, probe_positions, 'k-', 'LineWidth', 2);
xlabel('Sessions');
ylabel('Probe position (um)');

xlim(ax_motion, [0.5, n_session+0.5]);
xlabel(ax_session, 'Sessions');

%
similarity_all = [sum(Output.SimilarityAll.*Output.SimilarityWeights, 2), Output.SimilarityAll];
similarity_names = ['Weighted sum', Output.SimilarityNames];
similarity_matched = cell(1, length(similarity_names));
similarity_unmatched = cell(1, length(similarity_names));

is_matched = arrayfun(@(x)Output.ClusterMatrix(Output.SimilarityPairs(x,1), Output.SimilarityPairs(x,2)), 1:size(similarity_all,1));
idx_matched = find(is_matched == 1);
idx_unmatched = find(is_matched == 0);

for k = 1:length(similarity_names)
    similarity_matched{k} = similarity_all(idx_matched,k);
    similarity_unmatched{k} = similarity_all(idx_unmatched,k);
end

ax_hist = EasyPlot.createGridAxes(fig, 1, length(similarity_names),...
    'Width', 3,...
    'Height', 3,...
    'MarginLeft', 0.5);

EasyPlot.place(ax_hist, ax_session, 'right');
EasyPlot.align(ax_hist, ax_session, 'top');
EasyPlot.move(ax_hist, 'dx', 1.5);

for k = 1:length(ax_hist)
    histogram(ax_hist{k}, similarity_unmatched{k}, 'FaceColor', 'k', 'BinWidth', 0.1, 'Normalization', 'probability');
    histogram(ax_hist{k}, similarity_matched{k}, 'FaceColor', 'b', 'BinWidth', 0.1, 'Normalization', 'probability');
    xlabel(ax_hist{k}, similarity_names{k});

    if k > 1
        title(ax_hist{k}, ['weight = ', num2str(Output.SimilarityWeights(k-1), '%.3f')]);
    else
        xline(ax_hist{k}, Output.SimilarityThreshold, 'k:', 'lineWidth', 2);
    end
end

ylabel(ax_hist{1}, 'Probability');
EasyPlot.setYLim(ax_hist);
EasyPlot.setYLim(ax_hist, [0, ax_hist{1}.YLim(2)]);

% scatter
idx_scatter = cell(length(similarity_names) - 2, 1);
for k = 1:length(idx_scatter)
    idx_scatter{k} = [2, k+2];
end

ax_scatter = EasyPlot.createGridAxes(fig, 1, max(1, length(idx_scatter)),...
    'Width', 3,...
    'Height', 3,...
    'MarginLeft', 1,...
    'MarginRight', 1);

EasyPlot.align(ax_scatter, ax_hist, 'horizontalCenter');
EasyPlot.align(ax_scatter, ax_motion, 'top');

max_points = 5000;
for k = 1:length(idx_scatter)
    if length(similarity_unmatched{idx_scatter{k}(1)}) > max_points
        idx_plot = randperm(length(similarity_unmatched{idx_scatter{k}(1)}), max_points);
    else
        idx_plot = 1:length(similarity_unmatched{idx_scatter{k}(1)});
    end

    plot(ax_scatter{k},...
        similarity_unmatched{idx_scatter{k}(1)}(idx_plot),...
        similarity_unmatched{idx_scatter{k}(2)}(idx_plot), 'k.',...
        'MarkerSize', 1);

    if length(similarity_matched{idx_scatter{k}(1)}) > max_points
        idx_plot = randperm(length(similarity_matched{idx_scatter{k}(1)}), max_points);
    else
        idx_plot = 1:length(similarity_matched{idx_scatter{k}(1)});
    end

    plot(ax_scatter{k},...
        similarity_matched{idx_scatter{k}(1)}(idx_plot),...
        similarity_matched{idx_scatter{k}(2)}(idx_plot), 'b.',...
        'MarkerSize', 1);

    xlabel(ax_scatter{k}, similarity_names{idx_scatter{k}(1)});
    ylabel(ax_scatter{k}, similarity_names{idx_scatter{k}(2)});
end

EasyPlot.legend(ax_scatter{end}, {'Unmatched', 'Matched'},...
    'location', 'northeastoutside');

% plot the p(matched, delta_session)
n_matched_matrix = zeros(n_session);
for k = 1:n_cluster
    units = find(Output.IdxCluster == k);
    for j = 1:length(units)
        for i = j+1:length(units)
            n_matched_matrix(sessions(units(j)), sessions(units(i))) = n_matched_matrix(sessions(units(j)), sessions(units(i)))+1;
            n_matched_matrix(sessions(units(i)), sessions(units(j))) = n_matched_matrix(sessions(units(i)), sessions(units(j)))+1;
        end
    end
end

d_session = -n_session+1:n_session-1;
p_matched = cell(1, length(d_session));
p_matched_matrix = zeros(n_session);
for k = 1:n_session
    for j = k+1:n_session   
        p_matched_matrix(k,j) = n_matched_matrix(k,j)./n_units_each_session(k);
        p_matched_matrix(j,k) = n_matched_matrix(k,j)./n_units_each_session(j);
        
        idx_this = find(d_session == j-k);
        p_matched{idx_this} = [p_matched{idx_this}, n_matched_matrix(k,j)./n_units_each_session(j)];
        idx_this = find(d_session == k-j);
        p_matched{idx_this} = [p_matched{idx_this}, n_matched_matrix(k,j)./n_units_each_session(k)];
    end
end

p_matches_mean = zeros(1, length(d_session));
p_matches_std = zeros(1, length(d_session));
for k = 1:length(d_session)
    p_matches_mean(k) = mean(p_matched{k});
    p_matches_std(k) = std(p_matched{k});
end

ax_match_colormap = EasyPlot.createAxesAgainstAxes(fig, ax_motion, 'bottom',...
    'Width', 3,...
    'Height', 3,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

ax_p_matched = EasyPlot.createAxesAgainstAxes(fig, ax_match_colormap, 'right',...
    'Width', 3,...
    'Height', 3,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

EasyPlot.move(ax_p_matched, 'dx', 1.5);

imagesc(ax_match_colormap, p_matched_matrix);
xlim(ax_match_colormap, [0.5, n_session+0.5]);
ylim(ax_match_colormap, [0.5, n_session+0.5]);
xlabel(ax_match_colormap, 'Sessions');
ylabel(ax_match_colormap, 'Sessions');
EasyPlot.colorbar(ax_match_colormap, 'label', 'Prob. of matched');

plot(ax_p_matched, d_session, p_matches_mean, 'k.-');
% EasyPlot.plotShaded(ax_p_matched, 1:n_session-1, [p_matches_mean-p_matches_std; p_matches_mean+p_matches_std],...
%     'shadedColor', 'black');

xlabel(ax_p_matched, '\Delta session');
ylabel(ax_p_matched, 'Prob. of matched');
xlim(ax_p_matched, [-n_session+0.5, n_session-0.5]);

% Presence of units
[~, idx_sort_session] = sort(sessions*1e8 + locations(:,2)');
n_unique_units = n_cluster + sum(Output.IdxCluster == -1);
presence_matrix = zeros(n_unique_units, n_session);

cluster_count = 1;
cluster_id_plotted = [];
for k = 1:length(idx_sort_session)
    unit_this = idx_sort_session(k);
    if Output.IdxCluster(unit_this) == -1
        presence_matrix(cluster_count, sessions(unit_this)) = 1;
        cluster_count = cluster_count + 1;
        continue
    end

    cluster_id = Output.IdxCluster(unit_this);
    if any(cluster_id_plotted == cluster_id)
        continue
    end

    cluster_id_plotted = [cluster_id_plotted, cluster_id];
    presence_matrix(cluster_count, sessions(Output.IdxCluster == cluster_id)) = 1;
    cluster_count = cluster_count+1;
end

assert(size(presence_matrix, 1) == n_unique_units);

ax_presence = EasyPlot.createAxesAgainstAxes(fig, ax_match_colormap, 'bottom',...
    'Height', 10,...
    'MarginBottom', 1,...
    'MarginLeft', 1);
ax_presence.Position(3) = ax_p_matched.Position(1)+ax_p_matched.Position(3)-ax_match_colormap.Position(1);

imagesc(ax_presence, 1 - presence_matrix);
xlim(ax_presence, [0.5, n_session+0.5]);
ylim(ax_presence, [0.5, n_unique_units+0.5]);
EasyPlot.colormap(ax_presence, 'gray');
xlabel(ax_presence, 'Sessions');
ylabel(ax_presence, 'Units');

% result 
ax_colormap = EasyPlot.createGridAxes(fig, 2, 2,...
    'Width', 7,...
    'Height', 7,...
    'MarginTop', 0.2,...
    'MarginBottom', 0.2,...
    'MarginLeft', 0.2,...
    'MarginRight', 0.2);
EasyPlot.place(ax_colormap, ax_presence, 'right');
EasyPlot.align(ax_colormap, {ax_match_colormap, ax_p_matched, ax_presence}, 'verticalCenter');

EasyPlot.move(ax_colormap, 'dx', 1.5);

[~, idx_sort_session] = sort(sessions*1e8 + locations(:,2)');
n_unit_cumsum = cumsum(n_units_each_session);

imagesc(ax_colormap{1,1}, 1- Output.ClusterMatrix(idx_sort_session, idx_sort_session));
imagesc(ax_colormap{1,2}, Output.SimilarityMatrix(idx_sort_session, idx_sort_session));
imagesc(ax_colormap{2,1}, 1 - Output.ClusterMatrix(Output.IdxSort, Output.IdxSort));
imagesc(ax_colormap{2,2}, Output.SimilarityMatrix(Output.IdxSort, Output.IdxSort));

EasyPlot.colormap(ax_colormap(:,1), 'gray');

EasyPlot.setXLim(ax_colormap, [0.5, Output.NumUnits+0.5]);
EasyPlot.setYLim(ax_colormap, [0.5, Output.NumUnits+0.5]);

EasyPlot.setXTicksAndLabels(ax_colormap{1,1}, n_unit_cumsum, '');
EasyPlot.setYTicksAndLabels(ax_colormap{1,1}, n_unit_cumsum, [1, nan(1, n_session-2), n_session]);
EasyPlot.setXTicksAndLabels(ax_colormap{1,2}, n_unit_cumsum, '');
EasyPlot.setYTicksAndLabels(ax_colormap{1,2}, n_unit_cumsum, '');

x_ticks_units = ax_colormap{2,1}.XTick;
EasyPlot.setXTicksAndLabels(ax_colormap{2,1}, x_ticks_units, x_ticks_units);
EasyPlot.setYTicksAndLabels(ax_colormap{2,1}, x_ticks_units, '');
EasyPlot.setXTicksAndLabels(ax_colormap{2,2}, x_ticks_units, '');
EasyPlot.setYTicksAndLabels(ax_colormap{2,2}, x_ticks_units, '');

EasyPlot.ylabel(ax_colormap{1,1}, 'Sessions');
EasyPlot.xlabel(ax_colormap{2,1}, 'Units (sorted)');
title(ax_colormap{1,1}, 'Matches');
title(ax_colormap{1,2}, 'Similarity');

EasyPlot.colorbar(ax_colormap{end}, 'label', 'Similarity',...
    'MarginRight', 1);

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Overview'), 'dpi', 300);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Overview'), 'type', 'pdf', 'dpi', 300);

end