function waveforms_corrected = computeCorrectedWaveforms(user_settings, waveforms_all, channel_locations, sessions, locations, Motion)
% waveforms_all: n_unit x n_channel x n_sample array
% channel_locations: n_channel x 2 array
% sessions: n_unit x 1 array
% locations: n_unit x 2 array

n_unit = size(waveforms_all, 1);
n_channel = size(waveforms_all, 2);
n_sample = size(waveforms_all, 3);
n_templates = user_settings.waveformCorrection.n_templates;

waveforms_corrected = zeros(n_unit, n_channel, n_sample, n_templates);

min_channel_depth = min(channel_locations(:,2));
max_channel_depth = max(channel_locations(:,2));

motion_bottom = Motion.LinearScale*Motion.Linear*min_channel_depth + Motion.Constant;
motion_top = Motion.LinearScale*Motion.Linear*max_channel_depth + Motion.Constant;

min_motion = min([motion_bottom, motion_top]);
max_motion = max([motion_bottom, motion_top]);
fprintf('The range of motion: [%.1f μm ~ %.1f μm]\n', min_motion, max_motion);

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
    
    Motion_this = Motion;
    if n_templates == 2
        if i_template == 1
            Motion_this.Constant = Motion_this.Constant - min_motion;
        else
            Motion_this.Constant = Motion_this.Constant - max_motion;
        end
    end
    
    parfor k = 1:n_unit
        location_this = locations(k,:);

        dy = Motion_this.LinearScale*Motion_this.Linear(sessions(k))*location_this(2) + Motion_this.Constant(sessions(k));
        location_new = location_this;
        location_new(2) = location_new(2) - dy;
     
        waveforms_corrected(k,:,:,i_template) = waveformEstimation(...
            squeeze(waveforms_all(k,:,:)), location_this, channel_locations, location_new);
    
        updateParallel(1);
    end
    progBar.release();
end

end