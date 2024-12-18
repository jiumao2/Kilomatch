%% Recompute the similarities
max_distance = user_settings.clustering.max_distance;

similarity_waveform = zeros(1e7, 1);
similarity_raw_waveform = zeros(1e7, 1);
similarity_ISI = zeros(1e7, 1);
similarity_AutoCorr = zeros(1e7, 1);
similarity_PC = zeros(1e7, 1);
similarity_PETH = zeros(1e7, 1);
distance = zeros(1e7, 1);
idx_unit_pairs = zeros(1e7, 2);

count = 0;
% find the best unit pairs between days and the corresponding movement
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
        similarity_waveform(count) = waveformSimilarityMotionCorrected(waveforms_corrected([k,j],:,:), waveform_channels([k,j],:));
        similarity_raw_waveform(count) = waveformSimilarityMotionCorrected(waveforms([k,j],:,:), waveform_channels([k,j],:));
        similarity_ISI(count) = ISI_Similarity(spikeInfo(k), spikeInfo(j));
        similarity_AutoCorr(count) = autocorrelogramSimilarity(spikeInfo(k), spikeInfo(j));
        similarity_PC(count) = PC_SimilarityMotionCorrected(PC_features_corrected([k,j],:,:), PC_channels([k,j],:));
        similarity_PETH(count) = PETH_Similarity(spikeInfo(k), spikeInfo(j));
        distance_this = spikeInfo(j).Location - spikeInfo(k).Location;
        distance_this(2) = distance(2) - (positions(idx_block_j, session_j) - positions(idx_block_k, session_k));
        distance(count) = sqrt(sum(distance_this.^2));
        idx_unit_pairs(count,:) = [k,j];
    end
    if mod(k, 50) == 1
        toc
        fprintf('%d / %d done!\n', k, length(spikeInfo));
    end
end

similarity_waveform = similarity_waveform(1:count);
similarity_raw_waveform = similarity_raw_waveform(1:count);
similarity_ISI = similarity_ISI(1:count);
similarity_AutoCorr = similarity_AutoCorr(1:count);
similarity_PC = similarity_PC(1:count);
similarity_PETH = similarity_PETH(1:count);
distance = distance(1:count);

idx_unit_pairs = idx_unit_pairs(1:count,:);
session_pairs = [[spikeInfo(idx_unit_pairs(:,1)).SessionIndex]', [spikeInfo(idx_unit_pairs(:,2)).SessionIndex]'];

save(fullfile(user_settings.output_folder, 'AllSimilarity.mat'),...
    'similarity_waveform', 'similarity_raw_waveform', 'similarity_ISI', 'similarity_AutoCorr','similarity_PC', 'similarity_PETH',...
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

ax_PC = EasyPlot.createAxesAgainstAxes(fig, ax_autoCorr, 'right',...
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

histogram(ax_PC, similarity_PC, 'BinWidth', 0.2, 'Normalization', 'probability');
xlabel(ax_PC, 'PC similarity');
ylabel(ax_PC, 'Prob.');

EasyPlot.setYLim({ax_waveform, ax_autoCorr, ax_ISI, ax_PC});

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/AllSimilarity'));

