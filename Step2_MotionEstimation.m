% Get all features
if isfield(spikeInfo, 'ISI') &&...
        (any(strcmpi(user_settings.motionEstimation.features, 'ISI')) ||...
        any(strcmpi(user_settings.clustering.features, 'ISI')))
    ISI_features = cat(1, spikeInfo.ISI);
else
    ISI_features = [];
end

if isfield(spikeInfo, 'AutoCorr') &&...
        (any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr')) ||...
        any(strcmpi(user_settings.clustering.features, 'AutoCorr')))
    AutoCorr_features = cat(1, spikeInfo.AutoCorr);
else
    AutoCorr_features = [];
end

if isfield(spikeInfo, 'PETH') &&...
        (any(strcmpi(user_settings.motionEstimation.features, 'PETH')) ||...
        any(strcmpi(user_settings.clustering.features, 'PETH')))
    PETH_features = cat(1, spikeInfo.PETH);
else
    PETH_features = [];
end

%% Estimating the motion of the electrode
disp('---------------Motion Estimation---------------');
max_distance = user_settings.motionEstimation.max_distance;
n_nearest_channels = user_settings.waveformCorrection.n_nearest_channels;

unit_locations = zeros(1, length(spikeInfo));
for k = 1:length(unit_locations)
    unit_locations(k) = spikeInfo(k).Location(2);
end

y_distance_matrix = abs(unit_locations - unit_locations');
[idx_row, idx_col] = ind2sub(size(y_distance_matrix), 1:numel(y_distance_matrix));
idx_good = find(y_distance_matrix(:) <= max_distance & idx_col' > idx_row');
idx_unit_pairs = [idx_row(idx_good); idx_col(idx_good)]';

n_pairs = size(idx_unit_pairs, 1);

% compute waveform similarity
channel_locations = [spikeInfo(1).Xcoords, spikeInfo(1).Ycoords];
idx_nearest = knnsearch(channel_locations, channel_locations, 'K', n_nearest_channels);
idx_nearest_sorted = sort(idx_nearest, 2);

[idx_nearest_unique, ~, idx_groups] = unique(idx_nearest_sorted, 'rows');

waveform_all = zeros(length(spikeInfo), size(spikeInfo(1).Waveform, 1), size(spikeInfo(1).Waveform, 2));
waveform_similarity_matrix = zeros(length(spikeInfo));
for k = 1:length(spikeInfo)
    waveform_all(k,:,:) = spikeInfo(k).Waveform;
end

ptt = max(waveforms_corrected,[],3) - min(waveforms_corrected,[],3);
[~, ch] = max(ptt, [], 2);

progBar = ProgressBar(size(idx_nearest_unique, 1), ...
    'Title', 'Computing waveform features',...
    'UpdateRate', 1);

for k = 1:size(idx_nearest_unique, 1)
    idx_included = find(idx_groups == k);
    idx_units = find(arrayfun(@(x)any(idx_included==x), ch));

    if isempty(idx_units)
        progBar([], [], []);
        continue
    end

    waveform_this = reshape(waveform_all(:, idx_nearest_unique(k,:), :),...
        length(spikeInfo), []);

    temp = corrcoef(waveform_this');
    temp(isnan(temp)) = 0;
    temp = atanh(temp);
    
    waveform_similarity_matrix(idx_units,:) = temp(idx_units,:);

    progBar([], [], []);
end
progBar.release();

waveform_similarity_matrix = max(...
    cat(3, waveform_similarity_matrix, waveform_similarity_matrix'),...
    [], 3);

%% compute similarity
n_unit = length(spikeInfo);

similarity_waveform = arrayfun(...
        @(x)waveform_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';

similarity_ISI = zeros(n_pairs, 1);
ISI_similarity_matrix = zeros(n_unit);
if any(strcmpi(user_settings.motionEstimation.features, 'ISI'))
    ISI_similarity_matrix = corrcoef(ISI_features');
    ISI_similarity_matrix(isnan(ISI_similarity_matrix)) = 0;
    ISI_similarity_matrix = atanh(ISI_similarity_matrix);

    similarity_ISI = arrayfun(...
        @(x)ISI_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';
end

similarity_AutoCorr = zeros(n_pairs, 1);
AutoCorr_similarity_matrix = zeros(n_unit);
if any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr'))
    AutoCorr_similarity_matrix = corrcoef(AutoCorr_features');
    AutoCorr_similarity_matrix(isnan(AutoCorr_similarity_matrix)) = 0;
    AutoCorr_similarity_matrix = atanh(AutoCorr_similarity_matrix);

    similarity_AutoCorr = arrayfun(...
        @(x)AutoCorr_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';
end

similarity_PETH = zeros(n_pairs, 1);
PETH_similarity_matrix = zeros(n_unit);
if any(strcmpi(user_settings.motionEstimation.features, 'PETH'))
    PETH_similarity_matrix = corrcoef(PETH_features');
    PETH_similarity_matrix(isnan(PETH_similarity_matrix)) = 0;
    PETH_similarity_matrix = atanh(PETH_similarity_matrix);

    similarity_PETH = arrayfun(...
        @(x)PETH_similarity_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs)';
end

names_all = {'Waveform', 'ISI', 'AutoCorr', 'PETH'};
similarity_all = [similarity_waveform, similarity_ISI, similarity_AutoCorr, similarity_PETH];

similarity_matrix_all = cat(3,...
    waveform_similarity_matrix,...
    ISI_similarity_matrix,...
    AutoCorr_similarity_matrix,...
    PETH_similarity_matrix);

fprintf('Computing similarity done! Saved to %s ...\n', fullfile(user_settings.output_folder, 'SimilarityForCorretion.mat'));
toc;

% save the similarity
if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'SimilarityForCorretion.mat'),...
        'similarity_waveform', 'similarity_ISI', 'similarity_AutoCorr', 'similarity_PETH', 'idx_unit_pairs', '-nocompression');
end

% clear temp variables to save memory
clear temp ...
    similarity_waveform similarity_ISI similarity_AutoCorr similarity_PETH ...
    waveform_similarity_matrix ISI_similarity_matrix AutoCorr_similarity_matrix PETH_similarity_matrix...
    waveform_all unit_locations y_distance_matrix idx_row idx_col idx_good;

%% Pre-clustering
n_unit = length(spikeInfo);
sessions = [spikeInfo.SessionIndex];
n_session = max(sessions);

similarity_names = user_settings.motionEstimation.features';
idx_names = zeros(1, length(similarity_names));
for k = 1:length(similarity_names)
    idx_names(k) = find(strcmpi(names_all, similarity_names{k}));
end
similarity_all = similarity_all(:, idx_names);
similarity_matrix_all = similarity_matrix_all(:,:,idx_names);

weights = ones(1, length(similarity_names))./length(similarity_names);
similarity_matrix = squeeze(mean(similarity_matrix_all.*reshape(weights, 1, 1, n_features), 3));

% iterative clustering
for iter = 1:user_settings.motionEstimation.n_iter
    fprintf('Iteration %d starts!\n', iter);

    % HDBSCAN
    distance_matrix = 1./(1 + tanh(similarity_matrix));
    distance_matrix(eye(size(distance_matrix)) == 1) = 0;
    
    HDBSCAB_settings.min_samples = 1; % The number of samples in a neighborhood for a point to be considered as a core point.
    % This includes the point itself. When None, defaults to min_cluster_size.
    HDBSCAB_settings.cluster_selection_epsilon = 0; % A distance threshold. 
    % Clusters below this value will be merged. This is the minimum epsilon allowed.
    
    HDBSCAB_settings.min_cluster_size = 2;
    HDBSCAB_settings.max_cluster_size = n_session;
    HDBSCAB_settings.metric = 'precomputed';
    HDBSCAB_settings.data_folder = user_settings.output_folder;
    
    json_text = jsonencode(HDBSCAB_settings);
    
    fid = fopen(fullfile(user_settings.output_folder, 'HDBSCAN_settings.json'), 'w');
    fwrite(fid, json_text);
    fclose(fid);

    writeNPY(distance_matrix, fullfile(user_settings.output_folder, 'DistanceMatrix.npy'));
    
    % run HDBSCAN with Python
    system([fullfile(user_settings.path_to_python), ' ',...
        fullfile(path_kilomatch, 'Functions/main_hdbscan.py'), ' ',...
        fullfile(user_settings.output_folder, 'HDBSCAN_settings.json')]);
    
    idx_cluster_hdbscan = double(readNPY(fullfile(user_settings.output_folder, 'ClusterIndices.npy')));
    % MATLAB starts from 1
    idx_cluster_hdbscan(idx_cluster_hdbscan >= 0) = idx_cluster_hdbscan(idx_cluster_hdbscan >= 0)+1;
    
    n_cluster = max(idx_cluster_hdbscan);
    hdbscan_matrix = zeros(size(similarity_matrix), 'logical');
    for k = 1:n_cluster
        idx = find(idx_cluster_hdbscan == k);
        for j = 1:length(idx)
            for i = j+1:length(idx)
                hdbscan_matrix(idx(j), idx(i)) = 1;
                hdbscan_matrix(idx(i), idx(j)) = 1;
            end
        end
    end
    
    hdbscan_matrix(eye(n_unit) == 1) = 1;
    is_matched = arrayfun(@(x)hdbscan_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:n_pairs);

    if iter ~= user_settings.motionEstimation.n_iter        
        % LDA and update weights
        mdl = fitcdiscr(similarity_all, is_matched);
        temp = (mdl.Coeffs(1,2).Linear)';
        weights = temp./sum(temp);
        disp('Weights:');
        disp(strjoin(similarity_names, '   '));
        disp(weights);
        
        % update the similarity matrix
        similarity_matrix = squeeze(mean(similarity_matrix_all.*reshape(weights, 1, 1, n_features), 3));
    end
end

% set the threshold based on LDA results
similarity_thres = mdl.Coeffs(1,2).Const ./ (-mdl.Coeffs(1,2).Linear(1)) .* weights(1);

similarity = sum(similarity_all.*weights, 2);
idx_good = find(similarity' > similarity_thres & is_matched == 1);
n_pairs_included = length(idx_good);

% plot the final similarity score distribution
fig = EasyPlot.figure();
ax_similarity = EasyPlot.axes(fig,...
    'Width', 3,...
    'Height', 3,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

xline(ax_similarity, similarity_thres, 'k:', 'LineWidth', 2);

xlabel(ax_similarity, 'Waveform similarity');
ylabel(ax_similarity, 'Prob.');

histogram(ax_similarity, similarity, 'BinWidth', 0.2, 'Normalization', 'probability');
title(ax_similarity, [num2str(n_pairs_included), ' pairs are included']);
EasyPlot.cropFigure(fig);

if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/SimilarityThresholdForCorrection'));
end

fprintf('%d pairs of units are included for drift estimation!\n', n_pairs_included);

%% compute drift
nblock = user_settings.motionEstimation.n_block;
n_session = max([spikeInfo.SessionIndex]);
session_pairs = [[spikeInfo(idx_unit_pairs(:,1)).SessionIndex]', [spikeInfo(idx_unit_pairs(:,2)).SessionIndex]'];

% get all the good pairs and their distance
depth = zeros(1, length(idx_good));
dy = zeros(1, length(idx_good));
idx_1 = zeros(1, length(idx_good));
idx_2 = zeros(1, length(idx_good));
for k = 1:length(idx_good)
    unit1 = idx_unit_pairs(idx_good(k), 1);
    unit2 = idx_unit_pairs(idx_good(k), 2);
    d_this = mean([spikeInfo(unit2).Location(2), spikeInfo(unit1).Location(2)]);

    idx_1(k) = session_pairs(idx_good(k),1);
    idx_2(k) = session_pairs(idx_good(k),2);
    dy(k) = spikeInfo(unit2).Location(2) - spikeInfo(unit1).Location(2);
    depth(k) = d_this;
end

depth_edges = linspace(min(depth), max(depth), nblock+1);
depth_bins = 0.5*(depth_edges(1:end-1) + depth_edges(2:end));
idx_block = findNearestPoint(depth_bins, depth);

% compute the motion and 95CI
n_boot = 100;
positions = NaN(nblock, n_session);
positions_ci95 = zeros(2, nblock, n_session);
for k = 1:nblock
    dy_block = dy(idx_block == k);
    idx_1_block = idx_1(idx_block == k);
    idx_2_block = idx_2(idx_block == k);

    if length(unique([idx_1_block, idx_2_block])) ~= n_session
        disp('Some sessions are not included! Motion estimation failed!');
        continue
    end

    options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8);
    loss_fun = @(y) sum((dy_block - (y(idx_2_block) - y(idx_1_block))).^2, 'all');

    p = fminunc(loss_fun, rand(1, n_session), options);
    p = p - mean(p);
    positions(k,:) = p;

    options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8, 'Display', 'none');
    p_boot = zeros(n_boot, n_session);

    if isempty(gcp('nocreate'))
        parpool();
    end

    progBar = ProgressBar(n_boot, ...
        'IsParallel', true, ...
        'Title', 'Computing 95CI', ...
        'UpdateRate', 1 ...
        );
    progBar.setup([], [], []);

    parfor j = 1:n_boot
        idx_rand = randi(length(dy_block), 1, length(dy_block));
        dy_this = dy_block(idx_rand);
        idx_1_this = idx_1_block(idx_rand);
        idx_2_this = idx_2_block(idx_rand);
    
        loss_fun = @(y) sum((dy_this - (y(idx_2_this) - y(idx_1_this))).^2, 'all');
    
        p_this = fminunc(loss_fun, rand(1, n_session), options);
        p_boot(j,:) = p_this - mean(p_this);
    
        updateParallel(1);
    end
    progBar.release();

    p_ci95 = zeros(2, n_session);
    for j = 1:n_session
        p_ci95(1,j) = prctile(p_boot(:,j), 2.5);
        p_ci95(2,j) = prctile(p_boot(:,j), 97.5);
    end    

    positions_ci95(:,k,:) = p_ci95;
end

fprintf('%d / %d blocks are available!\n', sum(~isnan(positions(:,1))), nblock);

% interpolate the motion with nearest value if some blocks are not sampled
for k = 1:nblock
    if all(isnan(positions(k,:)))
        distance = abs(depth_bins-depth_bins(k));
        [~, idx_sort] = sort(distance);
        
        idx1 = NaN;
        d1 = NaN;
        for j = 1:length(idx_sort)
            if all(isnan(positions(idx_sort(j), :)))
                continue
            end

            if isnan(d1)
                d1 = distance(idx_sort(j));
                idx1 = idx_sort(j);
            elseif d1 == distance(idx_sort(j))
                positions(k,:) = mean(positions([idx1, idx_sort(j)], :));
            else
                positions(k,:) = positions(idx1, :);
            end
        end
    end
end

% Plot the drift
% compute the residues
dy_left = dy;
for k = 1:nblock
    dy_left(idx_block==k) = dy(idx_block==k) - (positions(k, idx_2(idx_block==k)) - positions(k, idx_1(idx_block==k)));
end

fig = EasyPlot.figure();
ax_drift = EasyPlot.axes(fig, 'Height', 3, 'Width', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 1);

ax_distance = EasyPlot.createAxesAgainstAxes(fig, ax_drift, 'right',...
    'Width', 8);

colors = [0.7,0.7,0.7; lines(1)];

if nblock > 1
    colors_block = parula(nblock);
else
    colors_block = 'k';
end

% drift
for k = 1:nblock
    plot(ax_drift, 1:n_session, positions(k,:), '-', 'Color', colors_block(k,:));
    EasyPlot.plotShaded(ax_drift, 1:n_session, squeeze(positions_ci95(:,k,:))', 'shadedColor', colors_block(k,:));
end

xlabel(ax_drift, 'Session');
ylabel(ax_drift, 'Motion (um)');

% raw distance
if nblock > 1
    for k = 1:nblock
        EasyPlot.plotShaded(ax_distance,...
            [depth_edges(k), depth_edges(k+1)],...
            [min(dy), max(dy); min(dy), max(dy)]', 'shadedColor', colors_block(k,:));
    end
end

plot(ax_distance, depth, dy, '.', 'Color', colors(1,:));
plot(ax_distance, depth, dy_left, '.', 'Color', colors(2,:));

xlabel(ax_distance, 'Depth (um)');
ylabel(ax_distance, 'Residues (um)');

xlim(ax_distance, [depth_edges(1), depth_edges(end)]);
ylim(ax_distance, [min(dy), max(dy)]);

EasyPlot.cropFigure(fig);
if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Motion'));
end

% save data
if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'Motion.mat'), 'positions', 'depth_bins', 'nblock', '-nocompression');
end

% clear temp variables
clear temp similarity similarity_all similarity_matrix distance_matrix hdbscan_matrix similarity_waveform similarity_ISI similarity_AutoCorr similarity_PETH;
