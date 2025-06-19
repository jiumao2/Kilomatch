function Output = saveToOutput(user_settings, spikeInfo,...
    idx_clusters, cluster_matrix, locations, leafOrder, ...
    similarity_matrix, similarity_all, idx_unit_pairs, similarity_names, weights, thres, good_matches_matrix,...
    sessions, Motion, idx_units,...
    curation_pairs, curation_types, curation_type_names, num_removal)

% get matched_pairs
[idx_row, idx_col] = ind2sub(size(cluster_matrix), find(cluster_matrix == 1));
idx_good = find(idx_col > idx_row);
matched_pairs = [idx_row(idx_good), idx_col(idx_good)];

Output = struct();
Output.NumClusters = max(idx_clusters);
Output.NumUnits = length(idx_clusters);
Output.IdxUnit = idx_units;
Output.Locations = locations; % NumUnits x 3;
Output.IdxSort = leafOrder;
Output.IdxCluster = idx_clusters;
Output.SimilarityMatrix = similarity_matrix; % the output similarity
Output.SimilarityAll = similarity_all;
Output.SimilarityPairs = idx_unit_pairs; % the pairs used to compute similarity_all
Output.SimilarityNames = similarity_names;
Output.SimilarityWeights = weights;
Output.SimilarityThreshold = thres;
Output.GoodMatchesMatrix = good_matches_matrix;
Output.ClusterMatrix = cluster_matrix;
Output.MatchedPairs = matched_pairs;

Output.CurationPairs = curation_pairs;
Output.CurationTypes = curation_types;
Output.CurationTypeNames = curation_type_names;
Output.CurationNumRemoval = num_removal;

Output.Params = user_settings;
Output.NumSession = max(sessions);
Output.Sessions = sessions;
Output.Motion = Motion;

if isfield(spikeInfo, 'Session')
    Output.SessionNames = {spikeInfo.Session};
end

Output.RunTime = toc;
Output.DateTime = datestr(datetime('now'));

save(fullfile(user_settings.output_folder, 'Output.mat'), 'Output', '-nocompression');
fprintf('Kilomatch done! Output is saved to %s!\n', fullfile(user_settings.output_folder, 'Output.mat'));
fprintf('Found %d clusters and %d matches from %d units during %d sessions!\n',...
    Output.NumClusters, size(Output.MatchedPairs, 1), Output.NumUnits, Output.NumSession);

end