function waveforms_corrected = computeCorrectedWaveforms(user_settings, waveforms_all, channel_locations, sessions, locations, motion)
% waveforms_all: n_unit x n_channel x n_sample array
% channel_locations: n_channel x 2 array
% sessions: n_unit x 1 array
% locations: n_unit x 2 array

n_unit = size(waveforms_all, 1);
n_channel = size(waveforms_all, 2);
n_sample = size(waveforms_all, 3);
n_templates = user_settings.waveformCorrection.n_templates;

waveforms_corrected = zeros(n_unit, n_channel, n_sample, n_templates);

% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end

for i_template = 1:n_templates
    progBar = ProgressBar(n_unit, ...
        'IsParallel', true, ...
        'Title', 'Computing waveform features', ...
        'UpdateRate', 1 ...
        );
    progBar.setup([], [], []);
    
    motion_this = motion;
    if n_templates == 2
        if i_template == 1
            motion_this = motion_this - min(motion_this);
        else
            motion_this = motion_this - max(motion_this);
        end
    end
    
    parfor k = 1:n_unit
        location_this = locations(k,:);
        dy = motion_this(sessions(k));
        location_new = location_this;
        location_new(2) = location_new(2) - dy;
     
        waveforms_corrected(k,:,:,i_template) = waveformEstimation(...
            squeeze(waveforms_all(k,:,:)), location_this, channel_locations, location_new);
    
        updateParallel(1);
    end
    progBar.release();
end

% Save the corrected waveforms
save(fullfile(user_settings.output_folder, 'Waveforms.mat'),...
    'waveforms_corrected', '-nocompression');

end