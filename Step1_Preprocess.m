%% load the data
load(user_settings.path_to_data);

%% validate the data
n_session = max([spikeInfo.SessionIndex]);
if n_session ~= length(unique([spikeInfo.SessionIndex]))
    error('SessionIndex should start from 1 and be coninuous without any gaps!');
end

%% preprocessing data
chanMap.xcoords = spikeInfo(1).Xcoords;
chanMap.ycoords = spikeInfo(1).Ycoords;
chanMap.kcoords = spikeInfo(1).Kcoords;
chanMap.connected = ones(1, length(spikeInfo(1).Xcoords));
for k = 1:length(spikeInfo)
    % compute the location of each unit
    [x, y, z, amp] = spikeLocation(spikeInfo(k).Waveform, chanMap,...
        user_settings.spikeLocation.n_nearest_channels,...
        user_settings.spikeLocation.location_algorithm);

    spikeInfo(k).Location = [x,y,z];
    spikeInfo(k).Amplitude = amp;
    
    % compute the peak channel
    [~, idx_max] = max(max(spikeInfo(k).Waveform, [], 2) - min(spikeInfo(k).Waveform, [], 2));
    spikeInfo(k).Channel = idx_max;
    
    % compute the autocorrelogram feauture
    spike_times = spikeInfo(k).SpikeTimes;
    window = user_settings.autoCorr.window; % ms   
    s = binTimings(spike_times, 1);
    
    [auto_corr, lag] = xcorr(s, s, window);
    auto_corr(lag==0)=0;
    
    spikeInfo(k).AutoCorr = auto_corr;

    % compute the ISI feature
    isi = diff(spike_times);
    isi_hist = histcounts(isi,...
        'BinLimits', [0, user_settings.ISI.window],...
        'BinWidth', user_settings.ISI.binwidth);
    isi_freq = isi_hist./sum(isi_hist);
    isi_freq_smoothed = smoothdata(isi_freq, 'gaussian', 5*user_settings.ISI.gaussian_sigma);
    spikeInfo(k).ISI = isi_freq_smoothed;


    if mod(k, 100) == 1
        toc;
        fprintf('%d / %d done!\n', k, length(spikeInfo));
    end
end

% make a folder to store the data
if ~exist(fullfile(user_settings.output_folder), 'dir')
    mkdir(fullfile(user_settings.output_folder));
end

% make a folder to store the figures
if ~exist(fullfile(user_settings.output_folder, 'Figures'), 'dir')
    mkdir(fullfile(user_settings.output_folder, 'Figures'));
end

% Save the preprocessed data
fprintf('Saving to %s...\n', fullfile(user_settings.output_folder, 'spikeInfo.mat'));
save(fullfile(user_settings.output_folder, 'spikeInfo.mat'), 'spikeInfo');

