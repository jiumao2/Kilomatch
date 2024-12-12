nPCs = user_settings.PC_similarity.n_PCs;
PCs = getPCs(spikeInfo, nPCs);

fig = EasyPlot.figure();
ax = EasyPlot.axes(fig,...
    'Width', 5,...
    'Height', 5,...
    'YAxisVisible', 'off',...
    'XAxisVisible', 'off');

plot(ax, 1:size(PCs, 2), PCs, '-', 'LineWidth', 2);
title(ax, 'Principle components');

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/PCs'));
%% Computing PC features with corrected waveforms
disp('Computing PC features with corrected waveforms...');
algorithm = user_settings.waveformCorrection.interpolate_algorithm;
n_nearest_channels_PC = user_settings.PC_similarity.n_channels_precomputed;

PC_features = zeros(length(spikeInfo), n_nearest_channels_PC, nPCs);
PC_features_corrected = zeros(length(spikeInfo), n_nearest_channels_PC, nPCs);
PC_channels = zeros(length(spikeInfo), n_nearest_channels_PC);

channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
chanMap.xcoords = spikeInfo(1).Xcoords;
chanMap.ycoords = spikeInfo(1).Ycoords;
chanMap.kcoords = spikeInfo(1).Kcoords;
chanMap.connected = ones(1, length(spikeInfo(1).Xcoords));

for k = 1:length(spikeInfo)
    % Do it with corrected waveforms
    location = spikeInfo(k).Location;
    idx_block = findNearestPoint(depth_bins, location(2));
    dy = positions(idx_block, spikeInfo(k).SessionIndex);
    location_new = location;
    location_new(2) = location_new(2) - dy;

    distance_to_location = sqrt(sum((channel_locations - location(1:2)).^2, 2));
    [~, idx_sort] = sort(distance_to_location);
    
    idx_included = idx_sort(1:n_nearest_channels_PC);
    PC_channels(k,:) = idx_included;

    for j = 1:n_nearest_channels_PC
        x = channel_locations(idx_included(j), 1);
        y = channel_locations(idx_included(j), 2);

        waveform_corrected = waveformEstimation(spikeInfo(k).Waveform, location, chanMap, location_new,...
            x, y, algorithm);
        for i = 1:nPCs
            PC_features(k,j,i) = sum(spikeInfo(k).Waveform(idx_included(j),:) .* PCs(i,:));
            PC_features_corrected(k,j,i) = sum(waveform_corrected .* PCs(i,:));
        end
    end
end

%% Save PC feautres
save(fullfile(user_settings.output_folder, 'PCs.mat'),...
    'PCs', 'PC_features', 'PC_features_corrected', 'PC_channels');

