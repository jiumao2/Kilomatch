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

%% preprocessing data
chanMap.xcoords = spikeInfo(1).Xcoords;
chanMap.ycoords = spikeInfo(1).Ycoords;
chanMap.kcoords = spikeInfo(1).Kcoords;
chanMap.connected = ones(1, length(spikeInfo(1).Xcoords));

n_unit = length(spikeInfo);
waveforms = zeros(n_unit, size(spikeInfo(1).Waveform, 1), size(spikeInfo(1).Waveform, 2));
for k = 1:length(spikeInfo)
    waveforms(k,:,:) = spikeInfo(k).Waveform;
end

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
    [x, y, z, amp] = spikeLocation(squeeze(waveforms(k,:,:)), chanMap,...
        user_settings.spikeLocation.n_nearest_channels,...
        user_settings.spikeLocation.location_algorithm);

    locations_all(k,:) = [x,y,z];
    amp_all(k) = amp;

    updateParallel(1);
end
progBar.release();

%% Compute spike times related features
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
    sigma = user_settings.autoCorr.gaussian_sigma;
    auto_corr_all = zeros(n_unit, 2*window+1);
    parfor k = 1:n_unit
        st_this = spike_times{k};
        st_this = st_this-st_this(1)+1;
        max_st = round(max(st_this))+1;
    
        s = zeros(1, max_st, 'logical');
        s(round(st_this)) = 1;
        
        [auto_corr, lag] = xcorr(s, s, window);
        auto_corr(lag==0)=0;
        auto_corr(lag>0) = smoothdata(auto_corr(lag>0), 'gaussian', 5*sigma);
        auto_corr(lag<0) = smoothdata(auto_corr(lag<0), 'gaussian', 5*sigma);
        auto_corr = auto_corr./max(auto_corr);
        
        auto_corr_all(k,:) = auto_corr;
    
        updateParallel(1);
    end
    progBar.release();
end
%%
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


%% Collect all preprocessed data
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

