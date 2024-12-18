% Overview of the clustering result
overviewResults(Output, user_settings);

% Visualization of the quality of each clusters
load(fullfile(user_settings.output_folder, 'spikeInfo.mat'));
waveformsAll = load(fullfile(user_settings.output_folder, 'Waveforms.mat'));
PC_all = load(fullfile(user_settings.output_folder, 'PCs.mat'));

for k = 1:Output.NumClusters
    visualizeCluster(Output, k, spikeInfo, waveformsAll, PC_all, user_settings);
    close all;
end






