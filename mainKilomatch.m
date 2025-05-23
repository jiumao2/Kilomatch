% Set the path to Kilomatch and settings
path_kilomatch = '.\Kilomatch'; % The path where Kilomatch is installed
path_settings = '.\settings.json'; % Please make sure the settings in the file are accurate

addpath(path_kilomatch);
addpath(genpath(fullfile(path_kilomatch, 'Functions')));

user_settings = jsonc.jsoncDecode(fileread(path_settings)); % Read the settings
tic;

%% Run Kilomatch
% load the data
fprintf('Loading %s...\n', user_settings.path_to_data);
load(user_settings.path_to_data);

spikeInfo = preprocessSpikeInfo(user_settings, spikeInfo);

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
Output = saveToOutput(user_settings, spikeInfo,...
    idx_cluster_hdbscan_curated, hdbscan_matrix_curated, locations, leafOrder, ...
    similarity_matrix, similarity_all, idx_unit_pairs, feature_names, weights, thres, good_matches_matrix,...
    sessions, motion);

% plot the result
overviewResults(user_settings, Output);
