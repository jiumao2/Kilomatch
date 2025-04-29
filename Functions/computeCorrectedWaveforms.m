function waveforms_corrected = computeCorrectedWaveforms(user_settings, waveforms_all, channel_locations, sessions, locations, motion)
% waveforms_all: n_unit x n_channel x n_sample array
% channel_locations: n_channel x 2 array
% sessions: n_unit x 1 array
% locations: n_unit x 2 array

n_unit = size(waveforms_all, 1);
n_channel = size(waveforms_all, 2);
n_sample = size(waveforms_all, 3);

waveforms_corrected = zeros(n_unit, n_channel, n_sample);

% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end
progBar = ProgressBar(n_unit, ...
    'IsParallel', true, ...
    'Title', 'Computing waveform features', ...
    'UpdateRate', 1 ...
    );
progBar.setup([], [], []);

parfor k = 1:n_unit
    location_this = locations(k,:);
    dy = motion(sessions(k));
    location_new = location_this;
    location_new(2) = location_new(2) - dy;
 
    waveforms_corrected(k,:,:) = waveformEstimation(...
        squeeze(waveforms_all(k,:,:)), location_this, channel_locations, location_new);

    updateParallel(1);
end
progBar.release();

% Save the corrected waveforms
save(fullfile(user_settings.output_folder, 'Waveforms.mat'),...
    'waveforms_corrected', '-nocompression');

end