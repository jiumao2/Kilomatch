function spikeInfo = preprocessSpikeInfo(user_settings, spikeInfo)
% spikeInfo: 1 x n_unit struct

% make a folder to store the data
if ~exist(fullfile(user_settings.output_folder), 'dir')
    mkdir(fullfile(user_settings.output_folder));
end
fprintf('The output will be saved to %s!\n', user_settings.output_folder);

% make a folder to store the figures
if ~exist(fullfile(user_settings.output_folder, 'Figures'), 'dir')
    mkdir(fullfile(user_settings.output_folder, 'Figures'));
end

% validate the data
n_session = max([spikeInfo.SessionIndex]);
if n_session ~= length(unique([spikeInfo.SessionIndex]))
    error('SessionIndex should start from 1 and be coninuous without any gaps!');
end

% preprocessing data
channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];

n_unit = length(spikeInfo);
waveforms_all = cat(3, spikeInfo.Waveform);
waveforms_all = permute(waveforms_all, [3,1,2]);

% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end
disp('Start preprocessing spikeInfo!');
progBar = ProgressBar(n_unit, ...
    'IsParallel', true, ...
    'Title', 'Compute unit locations', ...
    'UpdateRate', 1 ...
    );
progBar.setup([], [], []);

locations_all = zeros(n_unit, 3);
amp_all = zeros(n_unit, 1);
parfor k = 1:n_unit
    [x, y, z, amp] = spikeLocation(squeeze(waveforms_all(k,:,:)), channel_locations,...
        user_settings.spikeLocation.n_nearest_channels,...
        user_settings.spikeLocation.location_algorithm);

    locations_all(k,:) = [x,y,z];
    amp_all(k) = amp;

    updateParallel(1);
end
progBar.release();

% Compute spike times related features
spike_times = cell(1, n_unit);
for k = 1:n_unit
    spike_times{k} = sort(spikeInfo(k).SpikeTimes);
end

if any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr')) ||...
            any(strcmpi(user_settings.clustering.features, 'AutoCorr'))
    progBar = ProgressBar(n_unit, ...
        'IsParallel', true, ...
        'Title', 'Preprocessing AutoCorr', ...
        'UpdateRate', 1 ...
        );
    progBar.setup([], [], []);
    
    window = user_settings.autoCorr.window; % ms   
    binwidth = user_settings.autoCorr.binwidth; % ms   
    sigma = user_settings.autoCorr.gaussian_sigma;
    auto_corr_all = zeros(n_unit, 2*window+1);
    parfor k = 1:n_unit
        [auto_corr, lag] = computeAutoCorr(spike_times{k}, window, binwidth);

        auto_corr(lag>0) = smoothdata(auto_corr(lag>0), 'gaussian', 5*sigma);
        auto_corr(lag<0) = smoothdata(auto_corr(lag<0), 'gaussian', 5*sigma);
        auto_corr = auto_corr./max(auto_corr);
        
        auto_corr_all(k,:) = auto_corr;
    
        updateParallel(1);
    end
    progBar.release();
end

if any(strcmpi(user_settings.motionEstimation.features, 'ISI')) ||...
            any(strcmpi(user_settings.clustering.features, 'ISI'))
    progBar = ProgressBar(n_unit, ...
        'IsParallel', true, ...
        'Title', 'Preprocessing ISI', ...
        'UpdateRate', 1 ...
        );
    progBar.setup([], [], []);
    
    window = user_settings.ISI.window; % ms   
    sigma = user_settings.ISI.gaussian_sigma;
    binwidth = user_settings.ISI.binwidth;

    isi_all = zeros(n_unit, window);
    parfor k = 1:n_unit
        isi = diff(spike_times{k});
        isi_hist = histcounts(isi,...
            'BinLimits', [0, window],...
            'BinWidth', binwidth);
        isi_freq = isi_hist./sum(isi_hist);
        isi_freq_smoothed = smoothdata(isi_freq, 'gaussian', 5*sigma);
        isi_all(k,:) = isi_freq_smoothed;
    
        updateParallel(1);
    end
    progBar.release();
end

% Collect all preprocessed data
for k = 1:n_unit
    % get the location of each unit
    spikeInfo(k).Location = locations_all(k,:);
    spikeInfo(k).Amplitude = amp_all(k);
    
    % compute the peak channel
    [~, idx_max] = max(max(spikeInfo(k).Waveform, [], 2) - min(spikeInfo(k).Waveform, [], 2));
    spikeInfo(k).Channel = idx_max;

    % get the autocorrelogram feauture
    if any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr')) ||...
            any(strcmpi(user_settings.clustering.features, 'AutoCorr'))
        spikeInfo(k).AutoCorr = auto_corr_all(k,:);
    end

    % get the ISI feature
    if any(strcmpi(user_settings.motionEstimation.features, 'ISI')) ||...
            any(strcmpi(user_settings.clustering.features, 'ISI'))
        spikeInfo(k).ISI = isi_all(k,:);
    end
end

% Save the preprocessed data
fprintf('Saving to %s...\n', fullfile(user_settings.output_folder, 'spikeInfo.mat'));
save(fullfile(user_settings.output_folder, 'spikeInfo.mat'), 'spikeInfo', '-nocompression');

end






