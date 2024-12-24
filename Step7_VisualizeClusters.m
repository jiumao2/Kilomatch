% Overview of the clustering result
overviewResults(Output, user_settings);
close all;

% Visualization of the quality of each clusters
load(fullfile(user_settings.output_folder, 'spikeInfo.mat'));
waveformsAll = load(fullfile(user_settings.output_folder, 'Waveforms.mat'));

for k = 1:Output.NumClusters
    visualizeCluster(Output, k, spikeInfo, waveformsAll, user_settings);
    close all;
end





