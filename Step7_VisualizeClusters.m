% Overview of the clustering result
overviewResults(Output, user_settings);
close all;

% Visualization of the quality of each clusters
if user_settings.plot_clusters
    load(fullfile(user_settings.output_folder, 'Waveforms.mat'));
    for k = 1:Output.NumClusters
        visualizeCluster(Output, k, spikeInfo, waveforms_corrected, user_settings);
        close all;
    end
end




