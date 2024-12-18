[hdbscan_matrix_curated, idx_cluster_hdbscan_curated] = autoCuration(...
    hdbscan_matrix, idx_cluster_hdbscan, leafOrder, similarity_matrix, sessions,...
    user_settings);

matched_pairs_curated = [];
for k = 1:length(hdbscan_matrix_curated)
    for j = k+1:length(hdbscan_matrix_curated)
        if hdbscan_matrix_curated(k,j) == 1
            matched_pairs_curated = [matched_pairs_curated; k, j];
        end
    end
end

% Plot the final results
fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, 1, 2,...
    'Width', 12,...
    'Height', 12,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 0.2);

imagesc(ax_all{1}, hdbscan_matrix_curated(leafOrder,leafOrder));
imagesc(ax_all{2}, similarity_matrix(leafOrder,leafOrder));

EasyPlot.setCLim(ax_all{2}, [0, 4]);
h = EasyPlot.colorbar(ax_all{2},...
    'label', 'Similarity',...
    'MarginRight', 1);

EasyPlot.setXLim(ax_all, [0.5, length(leafOrder)+0.5]);
EasyPlot.setYLim(ax_all, [0.5, length(leafOrder)+0.5]);
title(ax_all{1}, 'Curated result');
title(ax_all{2}, 'Similarity matrix');

linkaxes([ax_all{1}, ax_all{2}]);

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/CuratedResult'));
savefig(fig, fullfile(user_settings.output_folder, 'Figures/CuratedResult.fig'));

%% save the curated result
save(fullfile(user_settings.output_folder, 'CurationResults.mat'),...
    'hdbscan_matrix_curated', 'idx_cluster_hdbscan_curated', 'matched_pairs_curated',...
    'similarity_matrix', 'sessions', 'leafOrder');


%% save the final output
% construct a similarity matrix for each feature
DistanceMatrix = NaN(length(idx_cluster_hdbscan_curated));
WaveformSimilarityMatrix = NaN(length(idx_cluster_hdbscan_curated));
RawWaveformSimilarityMatrix = NaN(length(idx_cluster_hdbscan_curated));
PC_SimilarityMatrix = NaN(length(idx_cluster_hdbscan_curated));
PETH_SimilarityMatrix = NaN(length(idx_cluster_hdbscan_curated));
ISI_SimilarityMatrix = NaN(length(idx_cluster_hdbscan_curated));
AutoCorrSimilalrityMatrix = NaN(length(idx_cluster_hdbscan_curated));

for k = 1:size(idx_unit_pairs, 1)
    DistanceMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = distance(k);
    DistanceMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = distance(k);

    WaveformSimilarityMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = similarity_waveform(k);
    WaveformSimilarityMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = similarity_waveform(k);

    RawWaveformSimilarityMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = similarity_raw_waveform(k);
    RawWaveformSimilarityMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = similarity_raw_waveform(k);

    PC_SimilarityMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = similarity_PC(k);
    PC_SimilarityMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = similarity_PC(k);

    PETH_SimilarityMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = similarity_PETH(k);
    PETH_SimilarityMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = similarity_PETH(k);

    ISI_SimilarityMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = similarity_ISI(k);
    ISI_SimilarityMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = similarity_ISI(k);

    AutoCorrSimilalrityMatrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = similarity_AutoCorr(k);
    AutoCorrSimilalrityMatrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = similarity_AutoCorr(k);
end

Output = struct();
Output.NumClusters = max(idx_cluster_hdbscan_curated);
Output.NumUnits = length(idx_cluster_hdbscan_curated);
Output.Locations = cat(1, spikeInfo.Location); % NumUnits x 3;
Output.IdxSort = leafOrder;
Output.IdxCluster = idx_cluster_hdbscan_curated;
Output.SimilarMatrix = similarity_matrix; % the output similarity
Output.ClusterMatrix = hdbscan_matrix_curated;
Output.MatchedPairs = matched_pairs_curated;
Output.Params = user_settings;
Output.NumSession = max(sessions);
Output.Sessions = sessions;
Output.Motion = positions;
Output.Nblock = nblock;
Output.PCs = PCs;

% The features used / unused in Kilomatch that might be useful in manual curation
Output.DistanceMatrix = DistanceMatrix;
Output.WaveformSimilarityMatrix = WaveformSimilarityMatrix;
Output.RawWaveformSimilarityMatrix = RawWaveformSimilarityMatrix;
Output.PC_SimilarityMatrix = PC_SimilarityMatrix;
Output.ISI_SimilarityMatrix = ISI_SimilarityMatrix;
Output.AutoCorrSimilalrityMatrix = AutoCorrSimilalrityMatrix;
Output.PETH_SimilarityMatrix = PETH_SimilarityMatrix;

save(fullfile(user_settings.output_folder, 'Output.mat'), 'Output');
fprintf('Kilomatch done! Output is saved to %s!\n', fullfile(user_settings.output_folder, 'Output.mat'));
fprintf('Found %d clusters and %d matches from %d units during %d sessions!\n',...
    Output.NumClusters, size(Output.MatchedPairs, 1), Output.NumUnits, Output.NumSession);

close all;

