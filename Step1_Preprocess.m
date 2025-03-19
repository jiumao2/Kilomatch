%% load the data
fprintf('Loading %s...\n', user_settings.path_to_data);
load(user_settings.path_to_data);

% make a folder to store the data
if ~exist(fullfile(user_settings.output_folder), 'dir')
    mkdir(fullfile(user_settings.output_folder));
end
fprintf('The output will be saved to %s!\n', user_settings.output_folder);

% make a folder to store the figures
if ~exist(fullfile(user_settings.output_folder, 'Figures'), 'dir')
    mkdir(fullfile(user_settings.output_folder, 'Figures'));
end

%% validate the data
n_session = max([spikeInfo.SessionIndex]);
if n_session ~= length(unique([spikeInfo.SessionIndex]))
    error('SessionIndex should start from 1 and be coninuous without any gaps!');
end

% check PETH
if ~isfield(spikeInfo(1), 'PETH')
    disp('PETH not found in SpikeInfo! Filled with zeros instead!');
    for k = 1:length(spikeInfo)
        spikeInfo(k).PETH = zeros(1, 5);
    end
end

%% preprocessing data
chanMap.xcoords = spikeInfo(1).Xcoords;
chanMap.ycoords = spikeInfo(1).Ycoords;
chanMap.kcoords = spikeInfo(1).Kcoords;
chanMap.connected = ones(1, length(spikeInfo(1).Xcoords));

% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end
disp('Start preprocessing spikeInfo!');
progBar = ProgressBar(length(spikeInfo), ...
    'IsParallel', true, ...
    'Title', 'Preprocessing spikeInfo', ...
    'UpdateRate', 1 ...
    );
progBar.setup([], [], []);

parfor k = 1:length(spikeInfo)
    % compute the location of each unit
    [x, y, z, amp] = spikeLocation(spikeInfo(k).Waveform, chanMap,...
        user_settings.spikeLocation.n_nearest_channels,...
        user_settings.spikeLocation.location_algorithm);

    spikeInfo(k).Location = [x,y,z];
    spikeInfo(k).Amplitude = amp;
    
    % compute the peak channel
    [~, idx_max] = max(max(spikeInfo(k).Waveform, [], 2) - min(spikeInfo(k).Waveform, [], 2));
    spikeInfo(k).Channel = idx_max;
    
    spike_times = spikeInfo(k).SpikeTimes;

    % compute the autocorrelogram feauture
    if any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr')) ||...
            any(strcmpi(user_settings.clustering.features, 'AutoCorr'))
        window = user_settings.autoCorr.window; % ms   
        s = binTimings(spike_times, 1);
        
        [auto_corr, lag] = xcorr(s, s, window);
        auto_corr(lag==0)=0;
        auto_corr(lag>0) = smoothdata(auto_corr(lag>0), 'gaussian', 5*user_settings.autoCorr.gaussian_sigma);
        auto_corr(lag<0) = smoothdata(auto_corr(lag<0), 'gaussian', 5*user_settings.autoCorr.gaussian_sigma);
        auto_corr = auto_corr./max(auto_corr);
        
        spikeInfo(k).AutoCorr = auto_corr;
    end

    % compute the ISI feature
    if any(strcmpi(user_settings.motionEstimation.features, 'ISI')) ||...
            any(strcmpi(user_settings.clustering.features, 'ISI'))
        isi = diff(spike_times);
        isi_hist = histcounts(isi,...
            'BinLimits', [0, user_settings.ISI.window],...
            'BinWidth', user_settings.ISI.binwidth);
        isi_freq = isi_hist./sum(isi_hist);
        isi_freq_smoothed = smoothdata(isi_freq, 'gaussian', 5*user_settings.ISI.gaussian_sigma);
        spikeInfo(k).ISI = isi_freq_smoothed;
    end

    updateParallel(1);
end
progBar.release();

% Save the preprocessed data
fprintf('Saving to %s...\n', fullfile(user_settings.output_folder, 'spikeInfo.mat'));
save(fullfile(user_settings.output_folder, 'spikeInfo.mat'), 'spikeInfo', '-nocompression');

