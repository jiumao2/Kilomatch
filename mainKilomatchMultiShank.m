% Set the path to Kilomatch and settings
path_kilomatch = '.\Kilomatch'; % The path where Kilomatch is installed
path_settings = '.\settings.json'; % Please make sure the settings are correct

addpath(path_kilomatch);
addpath(genpath(fullfile(path_kilomatch, 'Functions')));

user_settings = jsonc.jsoncDecode(fileread(path_settings)); % Read the settings
tic;

%% Run Kilomatch
Step1_Preprocess;

output_folder = user_settings.output_folder;
path_to_data = user_settings.path_to_data;

% Read the shank ID of each unit
shanks_data = arrayfun(@(x)x.Kcoords(x.Channel), spikeInfo);
shankIDs = unique(spikeInfo(1).Kcoords);

% Run Kilomatch in each shank individually
for i_shank = 1:length(shankIDs)
    % Reload the data and only consider the current shank
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

    % Run the remaining processes
    Step2_MotionEstimation;
    Step3_ComputeWaveformFeatures;
    Step4_ComputeAllSimilarity;
    Step5_IterativeClustering;
    Step6_AutoCuration;
    Step7_VisualizeClusters;
end

%% Combine the Output
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

