%% Computing the corrected waveforms
algorithm = user_settings.waveformCorrection.interpolate_algorithm;
n_channel = size(spikeInfo(1).Waveform, 1);
n_sample = size(spikeInfo(1).Waveform, 2);

waveforms_corrected = zeros(length(spikeInfo), n_channel, n_sample);

channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
chanMap.xcoords = spikeInfo(1).Xcoords;
chanMap.ycoords = spikeInfo(1).Ycoords;
chanMap.kcoords = spikeInfo(1).Kcoords;
chanMap.connected = ones(1, length(spikeInfo(1).Xcoords));

% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end
progBar = ProgressBar(length(spikeInfo), ...
    'IsParallel', true, ...
    'Title', 'Computing waveform features', ...
    'UpdateRate', 1 ...
    );
progBar.setup([], [], []);

parfor k = 1:length(spikeInfo)
    % Do it with corrected waveforms
    location = spikeInfo(k).Location;
    idx_block = findNearestPoint(depth_bins, location(2));
    dy = positions(idx_block, spikeInfo(k).SessionIndex);
    location_new = location;
    location_new(2) = location_new(2) - dy;

    for j = 1:n_channel
        x = channel_locations(j, 1);
        y = channel_locations(j, 2);
        
        waveforms_corrected(k,j,:) = waveformEstimation(...
            spikeInfo(k).Waveform, location, chanMap, location_new,...
            x, y, algorithm);
    end

    updateParallel(1);
end
progBar.release();

%% Save the corrected waveforms
save(fullfile(user_settings.output_folder, 'Waveforms.mat'),...
    'waveforms_corrected', '-nocompression');

