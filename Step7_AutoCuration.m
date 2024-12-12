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

%% save the curated result
save(fullfile(user_settings.output_folder, 'CurationResults.mat'),...
    'hdbscan_matrix_curated', 'idx_cluster_hdbscan_curated', 'matched_pairs_curated',...
    'similarity_matrix', 'sessions', 'merge_thres', 'reject_thres', 'leafOrder');


%% save the final output
Output = struct();
Output.NumClusters = max(idx_cluster_hdbscan_curated);
Output.NumUnits = length(idx_cluster_hdbscan_curated);
Output.IdxSort = leafOrder;
Output.IdxCluster = idx_cluster_hdbscan_curated;
Output.SimilarMatrix = similarity_matrix;
Output.ClusterMatrix = hdbscan_matrix_curated;
Output.MatchedPairs = matched_pairs_curated;
Output.Params = user_settings;
Output.NumSession = max(sessions);
Output.Sessions = sessions;
Output.Motion = positions;
Output.Nblock = nblock;
Output.PCs = PCs;

save(fullfile(user_settings.output_folder, 'Output.mat'), 'Output');
fprintf('Kilomatch done! Output is saved to %s!\n', fullfile(user_settings.output_folder, 'Output.mat'));
fprintf('Found %d clusters and %d matches from %d units during %d sessions!\n',...
    Output.NumClusters, size(Output.MatchedPairs, 1), Output.NumUnits, Output.NumSession);

