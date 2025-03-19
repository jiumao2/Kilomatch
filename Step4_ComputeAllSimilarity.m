%% Recompute the similarities
max_distance = user_settings.clustering.max_distance;

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

% clear temp variables to save memory
clear corrected_locations y_distance_matrix idx_row idx_col idx_good;

%%
similarity_waveform = zeros(n_pairs, 1);
similarity_raw_waveform = zeros(n_pairs, 1);
similarity_ISI = zeros(n_pairs, 1);
similarity_AutoCorr = zeros(n_pairs, 1);
similarity_PETH = zeros(n_pairs, 1);
distance = zeros(n_pairs, 1);

% compute similarity
% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end

disp('Start computing similarity!');
progBar = ProgressBar(n_pairs, ...
    'IsParallel', true, ...
    'Title', 'Computing Similarity', ...
    'UpdateRate', 1 ...
    );
progBar.setup([], [], []);

parfor k = 1:n_pairs
    idx_A = idx_unit_pairs(k,1);
    idx_B = idx_unit_pairs(k,2);

    session_A = session_pairs(k,1);
    session_B = session_pairs(k,2);

    idx_block_A = findNearestPoint(depth_bins, spikeInfo(idx_B).Location(2));
    idx_block_B = findNearestPoint(depth_bins, spikeInfo(idx_A).Location(2));
    
    if any(strcmpi(user_settings.clustering.features, 'Waveform'))
        similarity_waveform(k) = waveformSimilarityMotionCorrected(waveforms_corrected([idx_A,idx_B],:,:), waveform_channels([idx_A,idx_B],:),...
            user_settings.waveformCorrection.n_nearest_channels);
        similarity_raw_waveform(k) = waveformSimilarityMotionCorrected(waveforms([idx_A,idx_B],:,:), waveform_channels([idx_A,idx_B],:),...
            user_settings.waveformCorrection.n_nearest_channels);
    end
    
    if any(strcmpi(user_settings.clustering.features, 'ISI'))
        similarity_ISI(k) = ISI_Similarity(spikeInfo(idx_A), spikeInfo(idx_B));
    end

    if any(strcmpi(user_settings.clustering.features, 'AutoCorr'))
        similarity_AutoCorr(k) = autocorrelogramSimilarity(spikeInfo(idx_A), spikeInfo(idx_B));
    end

    if any(strcmpi(user_settings.clustering.features, 'PETH'))
        similarity_PETH(k) = PETH_Similarity(spikeInfo(idx_A), spikeInfo(idx_B));
    end
    
    distance_this = spikeInfo(idx_B).Location - spikeInfo(idx_A).Location;
    distance_this(2) = distance_this(2) - (positions(idx_block_B, session_B) - positions(idx_block_A, session_A));
    distance(k) = sqrt(sum(distance_this.^2));
    
    updateParallel(1);
end
progBar.release();

fprintf('Computing similarity done! Saved to %s ...\n', fullfile(user_settings.output_folder, 'AllSimilarity.mat'));
toc;

if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'AllSimilarity.mat'),...
        'similarity_waveform', 'similarity_raw_waveform', 'similarity_ISI', 'similarity_AutoCorr', 'similarity_PETH',...
        'distance', 'idx_unit_pairs', 'session_pairs', '-nocompression');
end

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

