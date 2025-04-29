% Set the path to Kilomatch and settings
path_kilomatch = '.\Kilomatch'; % The path where Kilomatch is installed
path_settings = '.\settings.json'; % Please make sure the settings in the file are accurate

addpath(path_kilomatch);
addpath(genpath(fullfile(path_kilomatch, 'Functions')));

user_settings = jsonc.jsoncDecode(fileread(path_settings)); % Read the settings
output_folder = user_settings.output_folder;
path_to_data = user_settings.path_to_data;
tic;

%% Run Kilomatch
% load the data
fprintf('Loading %s...\n', path_to_data);
load(path_to_data);

spikeInfo = preprocessSpikeInfo(user_settings, spikeInfo);

% read the shank ID of each unit
shanks_data = arrayfun(@(x)x.Kcoords(x.Channel), spikeInfo);
shankIDs = unique(spikeInfo(1).Kcoords);

% run Kilomatch in each shank individually
for i_shank = 1:length(shankIDs)
    % reload the data and only consider the current shank
    load(fullfile(output_folder, 'spikeInfo.mat'));
    shankID = shankIDs(i_shank);
    spikeInfo = spikeInfo(shanks_data == shankID);
    
    % set the new output folders
    user_settings.output_folder = fullfile(output_folder, ['Shank', num2str(shankID)]);
    if ~exist(fullfile(user_settings.output_folder), 'dir')
        mkdir(fullfile(user_settings.output_folder));
    end   
    if ~exist(fullfile(user_settings.output_folder, 'Figures'), 'dir')
        mkdir(fullfile(user_settings.output_folder, 'Figures'));
    end

    % get necessary information
    sessions = [spikeInfo.SessionIndex];
    channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
    locations = cat(1, spikeInfo.Location);
    [ISI_features, AutoCorr_features, PETH_features] = getAllFeatures(user_settings, spikeInfo);
    waveforms_all = cat(3, spikeInfo.Waveform);
    waveforms_all = permute(waveforms_all, [3,1,2]);
    
    % motion estimation
    features_all_motion_estimation = user_settings.motionEstimation.features;
    n_iter_motion_estimation = length(features_all_motion_estimation);
    motion = zeros(1, max(sessions)); % initialize the motion to zeros
    for i_iter = 1:n_iter_motion_estimation
        % find nearby pairs
        max_distance = user_settings.motionEstimation.max_distance;
        idx_unit_pairs = getNearbyPairs(max_distance, sessions, locations, motion);
        
        % compute similarity matrix 
        feature_names = features_all_motion_estimation{i_iter}';
        if i_iter == 1
            [similarity_matrix_all, feature_names_all] = computeAllSimilarityMatrix( ...
                user_settings, waveforms_all, channel_locations, ISI_features, AutoCorr_features, PETH_features, feature_names, idx_unit_pairs);
        else
            [similarity_matrix_all, feature_names_all] = computeAllSimilarityMatrix( ...
                user_settings, waveforms_corrected, channel_locations, ISI_features, AutoCorr_features, PETH_features, feature_names, idx_unit_pairs);
        end
    
        % iterative HDBSCAN
        idx_features = cellfun(@(x)find(strcmpi(feature_names_all, x)), feature_names);
        [hdbscan_matrix, ~, similarity_matrix, ~, ~, similarity_thres] = ...
            iterativeClustering(user_settings, path_kilomatch, similarity_matrix_all(:,:,idx_features), feature_names, idx_unit_pairs, sessions);
        
        % compute drift
        motion = computeMotion(user_settings, similarity_matrix, hdbscan_matrix, ...
            idx_unit_pairs, similarity_thres, sessions, locations);
    
        % compute corrected waveforms and save to Waveforms.mat
        waveforms_corrected = computeCorrectedWaveforms(user_settings, waveforms_all, channel_locations, sessions, locations, motion);
    end
    
    % final clustering
    % find nearby pairs
    max_distance = user_settings.clustering.max_distance;
    idx_unit_pairs = getNearbyPairs(max_distance, sessions, locations, motion);
    
    % compute similarity matrix
    [similarity_matrix_all, feature_names_all] = computeAllSimilarityMatrix( ...
        user_settings, waveforms_corrected, channel_locations, ISI_features, AutoCorr_features, PETH_features, feature_names, idx_unit_pairs);
    
    % iterative HDBSCAN
    feature_names = user_settings.clustering.features';
    idx_features = cellfun(@(x)find(strcmpi(feature_names_all, x)), feature_names);
    [hdbscan_matrix, idx_cluster_hdbscan, similarity_matrix, similarity_all, weights, thres, good_matches_matrix, leafOrder] = ...
        iterativeClustering(user_settings, path_kilomatch, similarity_matrix_all(:,:,idx_features), feature_names, idx_unit_pairs, sessions);
    
    % Auto-curate the result
    [hdbscan_matrix_curated, idx_cluster_hdbscan_curated] = autoCuration(...
        user_settings, hdbscan_matrix, idx_cluster_hdbscan, good_matches_matrix, ...
        sessions, similarity_matrix, leafOrder);
    
    % save the final output
    Output = saveToOutput(user_settings,...
        idx_cluster_hdbscan_curated, hdbscan_matrix_curated, locations, leafOrder, ...
        similarity_matrix, similarity_all, idx_unit_pairs, feature_names, weights, thres, good_matches_matrix,...
        sessions, motion);
    
    % plot the result
    overviewResults(user_settings, Output);
end

% Combine the Output
load(fullfile(output_folder, 'spikeInfo.mat'));
n_session = max([spikeInfo.SessionIndex]);

% Create an empty Output
Output = struct(...
    'NumCluster', 0,...
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
    'Params', user_settings,...
    'NumSession', n_session,...
    'Sessions', zeros(1, length(spikeInfo)),...
    'Motion', [],...
    'Nblock', user_settings.motionEstimation.n_block);

waveforms_corrected = zeros(length(spikeInfo), size(spikeInfo(1).Waveform, 1), size(spikeInfo(1).Waveform, 2));

n_cluster = 0;
n_units = 0;
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

    Output.Sessions(idx_units) = data.Output.Sessions;
    Output.Motion = [Output.Motion; data.Output.Motion];

    Output.RunTime = data.Output.RunTime; % save the run time of the final shank
    
    % update waveforms
    waveforms_corrected(idx_units,:,:) = data_waveforms.waveforms_corrected;

    n_cluster = n_cluster + data.Output.NumClusters;
end

Output.NumClusters = n_cluster;

% Save the combined output
fprintf('Saving Output to %s ...\n', fullfile(output_folder, 'Output.mat'));
save(fullfile(output_folder, 'Output.mat'), 'Output', '-nocompression');
save(fullfile(output_folder, 'Waveforms.mat'), 'waveforms_corrected', '-nocompression');



