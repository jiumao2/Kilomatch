function Output = mergeOutput(user_settings, spikeInfo, shanks_data, output_folder)
n_session = max([spikeInfo.SessionIndex]);
shankIDs = unique(shanks_data);

% Create an empty Output
Output = struct(...
    'IdxUnit', 1:length(spikeInfo),...
    'IdxShank', shanks_data,...
    'NumClusters', 0,...
    'NumUnits', length(spikeInfo),...
    'Locations', zeros(length(spikeInfo), 3),...
    'IdxSort', zeros(1, length(spikeInfo)),...
    'IdxCluster', zeros(length(spikeInfo), 1),...
    'SimilarityMatrix', zeros(length(spikeInfo)),...
    'SimilarityAll', [],...
    'SimilarityPairs', [],...
    'SimilarityNames', {user_settings.clustering.features'},...
    'SimilarityWeights', [],...
    'SimilarityThreshold', [],...
    'GoodMatchesMatrix', zeros(length(spikeInfo)),...
    'ClusterMatrix', zeros(length(spikeInfo)),...
    'MatchedPairs', [],...
    'CurationPairs', [],...
    'CurationTypes', [],...
    'CurationTypeNames', [],...
    'CurationNumRemoval', 0,...
    'Params', user_settings,...
    'NumSession', n_session,...
    'Sessions', zeros(1, length(spikeInfo)),...
    'SessionNames', [],...
    'Motion', struct('Linear', [], 'Constant', [], 'LinearScale', []),...
    'RunTime', [],...
    'DateTime', []);
Output.SessionNames = cell(1, length(spikeInfo));

waveforms_corrected = zeros(length(spikeInfo), size(spikeInfo(1).Waveform, 1), size(spikeInfo(1).Waveform, 2), user_settings.waveformCorrection.n_templates);

n_cluster = 0;
for i_shank = 1:length(shankIDs)
    shankID = shankIDs(i_shank);

    fprintf('Loading Output from shank %d ...\n', shankID);
    data = load(fullfile(output_folder, ['Shank', num2str(shankID)], 'Output.mat'));
    data_waveforms = load(fullfile(output_folder, ['Shank', num2str(shankID)], 'Waveforms.mat'));
    idx_units = find(shanks_data == shankID);

    Output.Locations(idx_units, :) = data.Output.Locations;
    Output.IdxSort(idx_units) = data.Output.IdxSort + n_cluster;

    Output.IdxCluster(idx_units) = data.Output.IdxCluster + n_cluster;
    Output.IdxCluster(idx_units(data.Output.IdxCluster == -1)) = -1;
    Output.SimilarityMatrix(idx_units, idx_units) = data.Output.SimilarityMatrix;
    Output.SimilarityAll = [Output.SimilarityAll; data.Output.SimilarityAll];

    similarity_pairs = arrayfun(@(x)idx_units(x), data.Output.SimilarityPairs);
    Output.SimilarityPairs = [Output.SimilarityPairs; similarity_pairs];

    Output.SimilarityWeights = [Output.SimilarityWeights; data.Output.SimilarityWeights];
    Output.SimilarityThreshold = [Output.SimilarityThreshold; data.Output.SimilarityThreshold];
    Output.GoodMatchesMatrix(idx_units, idx_units) = data.Output.GoodMatchesMatrix;
    Output.ClusterMatrix(idx_units, idx_units) = data.Output.ClusterMatrix;

    matched_pairs = arrayfun(@(x)idx_units(x), data.Output.MatchedPairs);
    Output.MatchedPairs = [Output.MatchedPairs; matched_pairs];

    curation_pairs = arrayfun(@(x)idx_units(x), data.Output.CurationPairs);
    Output.CurationPairs = [Output.CurationPairs; curation_pairs];
    Output.CurationTypes = [Output.CurationTypes, data.Output.CurationTypes];
    Output.CurationTypeNames = data.Output.CurationTypeNames;
    Output.CurationNumRemoval = Output.CurationNumRemoval + data.Output.CurationNumRemoval;

    Output.Sessions(idx_units) = data.Output.Sessions;
    if isfield(data.Output, 'SessionNames')
        Output.SessionNames(idx_units) = data.Output.SessionNames;
    end
    Output.Motion(i_shank) = data.Output.Motion;

    Output.RunTime = data.Output.RunTime; % save the run time of the final shank
    Output.DateTime = datestr(datetime('now'));
    
    % update waveforms
    waveforms_corrected(idx_units,:,:,:) = data_waveforms.waveforms_corrected;

    n_cluster = n_cluster + data.Output.NumClusters;
end

Output.NumClusters = n_cluster;

% Save the combined output
fprintf('Saving Output to %s ...\n', fullfile(output_folder, 'Output.mat'));
save(fullfile(output_folder, 'Output.mat'), 'Output', '-nocompression');
save(fullfile(output_folder, 'Waveforms.mat'), 'waveforms_corrected', '-nocompression');

end