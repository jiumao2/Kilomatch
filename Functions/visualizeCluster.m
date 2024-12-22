function visualizeCluster(Output, cluster_id, spikeInfo, waveformsAll, user_settings)
% (1) Positions on the probe vs. sessions
% (2) Overlapped waveforms
% (3) ISI, AutoCorr, PC, PETH
% (4) Positions on the 2D scatter

units = find(Output.IdxCluster == cluster_id);
sessions = Output.Sessions(units);
[~, idx_sort] = sort(sessions);
units = units(idx_sort);
sessions = sessions(idx_sort);

locations = Output.Locations(units, 1:2);
colors = winter(length(units));

waveforms = waveformsAll.waveforms(units,:,:);
waveform_channels = waveformsAll.waveform_channels(units,:);

ISI = cat(1, spikeInfo(units).ISI);
AutoCorr = cat(1, spikeInfo(units).AutoCorr);
AutoCorr(:, (size(AutoCorr, 2)+1)/2) = NaN;
PETH = cat(1, spikeInfo(units).PETH);

fig = EasyPlot.figure();

% depth
ax_depth = EasyPlot.axes(fig,...
    'Width', 5,...
    'Height', 3,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

plot(ax_depth, sessions, locations(:,2), 'k-');
for k = 1:length(sessions)
    plot(ax_depth, sessions(k), locations(k,2), '.', 'Color', colors(k,:),...
        'MarkerSize', 20);
end

title(ax_depth, ['n = ', num2str(length(units))]);
xlim(ax_depth, [min(sessions)-0.5, max(sessions)+0.5]);
ylim(ax_depth,...
    [min(mean(locations(:,2)) - user_settings.clustering.max_distance/2, min(locations(:,2) - 5)),...
    max(mean(locations(:,2)) + user_settings.clustering.max_distance/2, max(locations(:,2) + 5))]);
xlabel(ax_depth, 'Sessions');
ylabel(ax_depth, 'Depth (um)');

% waveforms
ch = mode(waveform_channels(:,1));
amplitude = max(max(squeeze(waveforms(:,1,:)), [], 2) - min(squeeze(waveforms(:,1,:)), [], 2));

n_channels = 20;
n_samples = size(waveforms, 3);
channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
distance_to_location = sqrt(sum((channel_locations - channel_locations(ch,:)).^2, 2));
[~, idx_sort] = sort(distance_to_location);
ch_included = idx_sort(1:n_channels);

waveforms_plot = zeros(length(units), n_channels, n_samples);
for k = 1:length(units)
    for j = 1:n_channels
        idx = waveform_channels(k,:) == ch_included(j);
        waveforms_plot(k,j,:) = waveforms(k, idx, :);
    end
end

x_plot = cell(1, length(units));
y_plot = cell(1, length(units));

samples_plot = (1:size(waveforms_plot, 3)) - (size(waveforms_plot, 3)+1)/2;
sample_scale = 1;
x_scale = 3;
y_scale = 1; % based on the amplitude
waveform_scale = 1./amplitude*50;
for j = 1:length(units)
    for k = 1:n_channels
        x = channel_locations(ch_included(k), 1);
        y = channel_locations(ch_included(k), 2);
    
        x_plot{j} = [x_plot{j}, x*x_scale+samples_plot*sample_scale, NaN];
        y_plot{j} = [y_plot{j}, y*y_scale+squeeze(waveforms_plot(j,k,:))' * waveform_scale, NaN];
    end
end

ax_waveform = EasyPlot.createAxesAgainstAxes(fig, ax_depth, 'bottom',...
    'Height', 8,...
    'YAxisVisible', 'off',...
    'XAxisVisible', 'off');

for k = 1:length(units)
    plot(ax_waveform, x_plot{k}, y_plot{k}, '-', 'Color', colors(k,:));
end 

EasyPlot.setYLim(ax_waveform, [min(y_plot{1}) - 10, max(y_plot{1}) + 10]);
EasyPlot.setXLim(ax_waveform, [min(x_plot{1}) - 10, max(x_plot{1}) + 10]);

h_scalebar = EasyPlot.scalebar(ax_waveform, 'XY',...
    'location', 'southeast',...
    'xBarLabel', '1 ms',...
    'xBarRatio', 1,...
    'xBarLength', 30,...
    'yBarLabel', '0.2 mV',...
    'yBarRatio', waveform_scale,...
    'yBarLength', 100);

EasyPlot.move(h_scalebar, 'dx', 1);

% similarity
overlaps = {...
    [],...
    squeeze(waveforms_plot(:,1,:)),...
    ISI,...
    AutoCorr,...
    PETH};
similarity_matrix = {...
    Output.SimilarityMatrix(units, units),...
    Output.WaveformSimilarityMatrix(units, units),...
    Output.ISI_SimilarityMatrix(units, units),...
    Output.AutoCorrSimilalrityMatrix(units, units),...
    Output.PETH_SimilarityMatrix(units, units),...
    };

similarity_names = ['Weighted sum', Output.SimilarityNames];

ax_similarity = EasyPlot.createGridAxes(fig, 2, length(similarity_names),...
    'Width', 3,...
    'Height', 3,...
    'MarginLeft', 1,...
    'MarginBottom', 1);

EasyPlot.place(ax_similarity, ax_depth, 'right');
EasyPlot.align(ax_similarity, ax_depth, 'top');


EasyPlot.set(ax_similarity(2,:),...
    'Color', [0.7,0.7,0.7]);

for k = 1:length(ax_similarity)
    similarity_matrix{k}(eye(size(similarity_matrix{k})) == 1) = NaN;
    
    for j = 1:size(overlaps{k}, 1)
        if k > 1
            plot(ax_similarity{1,k}, 1:size(overlaps{k},2), overlaps{k}(j,:), '-', 'Color', colors(j,:));
        end
    end

    imagesc(ax_similarity{2,k}, similarity_matrix{k}, 'AlphaData', ~isnan(similarity_matrix{k}));
    title(ax_similarity{1,k}, similarity_names{k});

    if k == 1
        title(ax_similarity{2,k}, similarity_names{k});
    else
        title(ax_similarity{2,k}, ['weight = ', num2str(Output.SimilarityWeights(k-1), '%.3f')]);
    end
end

EasyPlot.setXLim(ax_similarity(2,:), [0.5, length(units)+0.5]);
EasyPlot.setYLim(ax_similarity(2,:), [0.5, length(units)+0.5]);
EasyPlot.setCLim(ax_similarity(2,:), [0, 4]);
EasyPlot.colorbar(ax_similarity{2,end}, 'label', 'Similarity', 'MarginRight', 1);

EasyPlot.xticklabels(ax_similarity(2,:), '');
EasyPlot.yticklabels(ax_similarity(2,:), '');
xlabel(ax_similarity{2,1}, 'Units (sorted by sessions)');

% 2D scatter
similarity_all = Output.SimilarityAll;
similarity_names = {'Waveform', 'ISI', 'Auto correlogram', 'PETH'};
similarity_matched = cell(1, length(similarity_names));
similarity_unmatched = cell(1, length(similarity_names));
similarity_this = cell(1, length(similarity_names));

is_matched = arrayfun(@(x)Output.ClusterMatrix(Output.SimilarityPairs(x,1), Output.SimilarityPairs(x,2)), 1:size(similarity_all,1));
idx_matched = find(is_matched == 1);
idx_unmatched = find(is_matched == 0);

temp_matrix = zeros(Output.NumUnits);
temp_matrix(units, units) = 1;
is_this_cluster = arrayfun(@(x)temp_matrix(Output.SimilarityPairs(x,1), Output.SimilarityPairs(x,2)), 1:size(similarity_all,1));
idx_this_cluster = find(is_this_cluster == 1);

similarity = sum(Output.SimilarityAll.*Output.SimilarityWeights, 2);
histogram(ax_similarity{1,1}, similarity(idx_unmatched), 'FaceColor', 'k', 'BinWidth', 0.2, 'Normalization', 'probability');
histogram(ax_similarity{1,1}, similarity(idx_matched), 'FaceColor', 'b', 'BinWidth', 0.2, 'Normalization', 'probability');
plot(ax_similarity{1,1}, similarity(idx_this_cluster), zeros(length(idx_this_cluster), 1), 'g.', 'MarkerSize', 3);

for k = 1:length(similarity_names)
    similarity_matched{k} = similarity_all(idx_matched,k);
    similarity_unmatched{k} = similarity_all(idx_unmatched,k);
    similarity_this{k} = similarity_all(idx_this_cluster,k);
end

idx_scatter = cell(1, length(similarity_names) - 1);
for k = 1:length(idx_scatter)
    idx_scatter{k} = [1, k+1];
end

ax_scatter = EasyPlot.createGridAxes(fig, 1, length(idx_scatter),...
    'Width', 4,...
    'Height', 4,...
    'MarginLeft', 1,...
    'MarginRight', 1,...
    'MarginBottom', 1);

EasyPlot.align(ax_scatter, ax_similarity, 'horizontalCenter');
EasyPlot.align(ax_scatter, ax_waveform, 'bottom');

max_points = 5000;
for k = 1:length(ax_scatter)
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

    plot(ax_scatter{k},...
        similarity_this{idx_scatter{k}(1)},...
        similarity_this{idx_scatter{k}(2)}, 'g+',...
        'MarkerSize', 5);

    xlabel(ax_scatter{k}, similarity_names{idx_scatter{k}(1)});
    ylabel(ax_scatter{k}, similarity_names{idx_scatter{k}(2)});
end

EasyPlot.legend(ax_scatter{end}, {'Matched', 'Unmatched', 'This cluster'},...
    'location', 'northeastoutside');

h = EasyPlot.setGeneralTitle([{ax_depth}, ax_similarity(1,:)], ['Cluster #', num2str(cluster_id)],...
    'fontSize', 15);
EasyPlot.move(h, 'dy', 1);

EasyPlot.cropFigure(fig);

if ~exist(fullfile(user_settings.output_folder, 'Figures/Clusters/'), 'dir')
    mkdir(fullfile(user_settings.output_folder, 'Figures/Clusters/'));
end
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Clusters', ['Cluster', num2str(cluster_id)]));

end