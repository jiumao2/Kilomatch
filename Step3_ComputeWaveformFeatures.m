%% Computing the corrected waveforms
algorithm = user_settings.waveformCorrection.interpolate_algorithm;
n_nearest_channels = user_settings.waveformCorrection.n_channels_precomputed;
n_sample = size(spikeInfo(1).Waveform, 2);

waveforms = zeros(length(spikeInfo), n_nearest_channels, n_sample);
waveforms_corrected = zeros(length(spikeInfo), n_nearest_channels, n_sample);
waveform_channels = zeros(length(spikeInfo), n_nearest_channels);

channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
chanMap.xcoords = spikeInfo(1).Xcoords;
chanMap.ycoords = spikeInfo(1).Ycoords;
chanMap.kcoords = spikeInfo(1).Kcoords;
chanMap.connected = ones(1, length(spikeInfo(1).Xcoords));

progBar = ProgressBar(...
    length(spikeInfo), ...
    'Title', 'Computing corrected waveforms' ...
    );
for k = 1:length(spikeInfo)
    % Do it with corrected waveforms
    location = spikeInfo(k).Location;
    idx_block = findNearestPoint(depth_bins, location(2));
    dy = positions(idx_block, spikeInfo(k).SessionIndex);
    location_new = location;
    location_new(2) = location_new(2) - dy;

    distance_to_location = sqrt(sum((channel_locations - location_new(1:2)).^2, 2));
    [~, idx_sort] = sort(distance_to_location);
    
    idx_included = idx_sort(1:n_nearest_channels);
    waveform_channels(k,:) = idx_included;

    for j = 1:n_nearest_channels
        x = channel_locations(idx_included(j), 1);
        y = channel_locations(idx_included(j), 2);
        
        waveforms(k,j,:) = spikeInfo(k).Waveform(idx_included(j),:);
        waveforms_corrected(k,j,:) = waveformEstimation(spikeInfo(k).Waveform, location, chanMap, location_new,...
            x, y, algorithm);
    end

    progBar([], [], []);
end
progBar.release();

%% Save the corrected waveforms
save(fullfile(user_settings.output_folder, 'Waveforms.mat'),...
    'waveforms', 'waveforms_corrected', 'waveform_channels', '-nocompression');

