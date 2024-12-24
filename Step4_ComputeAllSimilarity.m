%% Recompute the similarities
max_distance = user_settings.clustering.max_distance;

idx_unit_pairs = zeros(1e8, 2);
count = 0;

progBar = ProgressBar(...
    length(spikeInfo), ...
    'Title', 'Getting unit pairs' ...
    );
for k = 1:length(spikeInfo)
    for j = k+1:length(spikeInfo)
        session_k = spikeInfo(k).SessionIndex;
        session_j = spikeInfo(j).SessionIndex;

        idx_block_j = findNearestPoint(depth_bins, spikeInfo(j).Location(2));
        idx_block_k = findNearestPoint(depth_bins, spikeInfo(k).Location(2));
        
        if abs(spikeInfo(j).Location(2) - spikeInfo(k).Location(2) -...
                (positions(idx_block_j, session_j) - positions(idx_block_k, session_k))) > max_distance
            continue
        end
        
        count = count+1;
        idx_unit_pairs(count,:) = [k,j];
    end

    progBar([], [], []);
end
progBar.release();

idx_unit_pairs = idx_unit_pairs(1:count,:);
session_pairs = [[spikeInfo(idx_unit_pairs(:,1)).SessionIndex]', [spikeInfo(idx_unit_pairs(:,2)).SessionIndex]'];
n_pairs = size(idx_unit_pairs, 1);

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
    
    similarity_waveform(k) = waveformSimilarityMotionCorrected(waveforms_corrected([idx_A,idx_B],:,:), waveform_channels([idx_A,idx_B],:));
    similarity_raw_waveform(k) = waveformSimilarityMotionCorrected(waveforms([idx_A,idx_B],:,:), waveform_channels([idx_A,idx_B],:));
    similarity_ISI(k) = ISI_Similarity(spikeInfo(idx_A), spikeInfo(idx_B));
    similarity_AutoCorr(k) = autocorrelogramSimilarity(spikeInfo(idx_A), spikeInfo(idx_B));
    similarity_PETH(k) = PETH_Similarity(spikeInfo(idx_A), spikeInfo(idx_B));
    distance_this = spikeInfo(idx_B).Location - spikeInfo(idx_A).Location;
    distance_this(2) = distance_this(2) - (positions(idx_block_B, session_B) - positions(idx_block_A, session_A));
    distance(k) = sqrt(sum(distance_this.^2));
    
    updateParallel(1);
end
progBar.release();

fprintf('Computing similarity done! Saved to %s ...\n', fullfile(user_settings.output_folder, 'AllSimilarity.mat'));
toc;

save(fullfile(user_settings.output_folder, 'AllSimilarity.mat'),...
    'similarity_waveform', 'similarity_raw_waveform', 'similarity_ISI', 'similarity_AutoCorr', 'similarity_PETH',...
    'distance', 'idx_unit_pairs', 'session_pairs');

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
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/AllSimilarity'));

