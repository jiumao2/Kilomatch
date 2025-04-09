%% Recompute the similarities
max_distance = user_settings.clustering.max_distance;
n_nearest_channels = user_settings.waveformCorrection.n_nearest_channels;

corrected_locations = zeros(1, length(spikeInfo));
for k = 1:length(corrected_locations)
    idx_block = findNearestPoint(depth_bins, spikeInfo(k).Location(2));
    corrected_locations(k) = spikeInfo(k).Location(2) - positions(idx_block, spikeInfo(k).SessionIndex);
end

y_distance_matrix = abs(corrected_locations - corrected_locations');

idx_col = floor((0:numel(y_distance_matrix)-1) ./ size(y_distance_matrix, 1))' + 1;
idx_row = mod((0:numel(y_distance_matrix)-1), size(y_distance_matrix, 1))' + 1;
idx_good = find(y_distance_matrix(:) <= max_distance & idx_col > idx_row);
idx_unit_pairs = [idx_row(idx_good), idx_col(idx_good)];

session_pairs = [[spikeInfo(idx_unit_pairs(:,1)).SessionIndex]', [spikeInfo(idx_unit_pairs(:,2)).SessionIndex]'];
n_pairs = size(idx_unit_pairs, 1);

% compute waveform similarity
channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
idx_nearest = knnsearch(channel_locations, channel_locations, 'K', n_nearest_channels);
idx_nearest_sorted = sort(idx_nearest, 2);

[idx_nearest_unique, ~, idx_groups] = unique(idx_nearest_sorted, 'rows');

waveform_similarity_matrix_corrected = zeros(length(spikeInfo));

ptt = max(waveforms_corrected,[],3) - min(waveforms_corrected,[],3);
[~, ch] = max(ptt, [], 2);

progBar = ProgressBar(size(idx_nearest_unique, 1), ...
    'Title', 'Computing waveform features', ...
    'UpdateRate', 1);
for k = 1:size(idx_nearest_unique, 1)
    idx_included = find(idx_groups == k);
    idx_units = find(arrayfun(@(x)any(idx_included==x), ch));

    if isempty(idx_units)
        progBar([], [], []);
        continue
    end

    waveform_this = reshape(waveforms_corrected(:, idx_nearest_unique(k,:), :),...
        length(spikeInfo), []);

    temp = corrcoef(waveform_this');
    temp(isnan(temp)) = 0;
    temp = atanh(temp);
    
    waveform_similarity_matrix_corrected(idx_units,:) = temp(idx_units,:);

    progBar([], [], []);
end
progBar.release();

waveform_similarity_matrix_corrected = max(...
    cat(3, waveform_similarity_matrix_corrected, waveform_similarity_matrix_corrected'),...
    [], 3);

n_unit = length(spikeInfo);

similarity_waveform = arrayfun(...
        @(x)waveform_similarity_matrix_corrected(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';

similarity_ISI = zeros(n_pairs, 1);
ISI_similarity_matrix = zeros(n_unit);
if any(strcmpi(user_settings.motionEstimation.features, 'ISI'))
    ISI_similarity_matrix = corrcoef(ISI_features');
    ISI_similarity_matrix(isnan(ISI_similarity_matrix)) = 0;
    ISI_similarity_matrix = atanh(ISI_similarity_matrix);

    similarity_ISI = arrayfun(...
        @(x)ISI_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';
end

similarity_AutoCorr = zeros(n_pairs, 1);
AutoCorr_similarity_matrix = zeros(n_unit);
if any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr'))
    AutoCorr_similarity_matrix = corrcoef(AutoCorr_features');
    AutoCorr_similarity_matrix(isnan(AutoCorr_similarity_matrix)) = 0;
    AutoCorr_similarity_matrix = atanh(AutoCorr_similarity_matrix);

    similarity_AutoCorr = arrayfun(...
        @(x)AutoCorr_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';
end

similarity_PETH = zeros(n_pairs, 1);
PETH_similarity_matrix = zeros(n_unit);
if any(strcmpi(user_settings.motionEstimation.features, 'PETH'))
    PETH_similarity_matrix = corrcoef(PETH_features');
    PETH_similarity_matrix(isnan(PETH_similarity_matrix)) = 0;
    PETH_similarity_matrix = atanh(PETH_similarity_matrix);

    similarity_PETH = arrayfun(...
        @(x)PETH_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';
end

similarity_matrix_all = cat(3,...
    waveform_similarity_matrix_corrected,...
    ISI_similarity_matrix,...
    AutoCorr_similarity_matrix,...
    PETH_similarity_matrix);

fprintf('Computing similarity done! Saved to %s ...\n', fullfile(user_settings.output_folder, 'AllSimilarity.mat'));
toc;

if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'AllSimilarity.mat'),...
        'similarity_waveform', 'similarity_ISI', 'similarity_AutoCorr', 'similarity_PETH',...
        'idx_unit_pairs', 'session_pairs', '-nocompression');
end

% clear temp variables to save memory
clear temp waveform_similarity_matrix_corrected ISI_similarity_matrix AutoCorr_similarity_matrix PETH_similarity_matrix...
    waveforms_corrected corrected_locations y_distance_matrix idx_row idx_col idx_good;
%%
fig = EasyPlot.figure();
ax_waveform = EasyPlot.axes(fig,...
    'Width', 5,...
    'Height', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

ax_ISI = EasyPlot.createAxesAgainstAxes(fig, ax_waveform, 'right',...
    'Width', 5,...
    'Height', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

ax_autoCorr = EasyPlot.createAxesAgainstAxes(fig, ax_ISI, 'right',...
    'Width', 5,...
    'Height', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

ax_PETH = EasyPlot.createAxesAgainstAxes(fig, ax_autoCorr, 'right',...
    'Width', 5,...
    'Height', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

histogram(ax_waveform, similarity_waveform, 'BinWidth', 0.2, 'Normalization', 'probability');
xlabel(ax_waveform, 'Waveform similarity');
ylabel(ax_waveform, 'Prob.');

histogram(ax_ISI, similarity_ISI, 'BinWidth', 0.2, 'Normalization', 'probability');
xlabel(ax_ISI, 'ISI similarity');
ylabel(ax_ISI, 'Prob.');

histogram(ax_autoCorr, similarity_AutoCorr, 'BinWidth', 0.2, 'Normalization', 'probability');
xlabel(ax_autoCorr, 'AutoCorr similarity');
ylabel(ax_autoCorr, 'Prob.');

histogram(ax_PETH, similarity_PETH, 'BinWidth', 0.2, 'Normalization', 'probability');
xlabel(ax_PETH, 'PETH similarity');
ylabel(ax_PETH, 'Prob.');

EasyPlot.setYLim({ax_waveform, ax_autoCorr, ax_ISI, ax_PETH});

EasyPlot.cropFigure(fig);

if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/AllSimilarity'));
end

