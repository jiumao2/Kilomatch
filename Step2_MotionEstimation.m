% Estimating the motion of the electrode
disp('---------------Motion Estimation---------------');
max_distance = user_settings.motionEstimation.max_motion_distance;

unit_locations = zeros(1, length(spikeInfo));
for k = 1:length(unit_locations)
    unit_locations(k) = spikeInfo(k).Location(2);
end

y_distance_matrix = abs(unit_locations - unit_locations');

idx_col = floor((0:numel(y_distance_matrix)-1) ./ size(y_distance_matrix, 1))' + 1;
idx_row = mod((0:numel(y_distance_matrix)-1), size(y_distance_matrix, 1))' + 1;
idx_good = find(y_distance_matrix(:) <= max_distance & idx_col > idx_row);
idx_unit_pairs = [idx_row(idx_good), idx_col(idx_good)];

n_pairs = size(idx_unit_pairs, 1);

% clear temp variables to save memory
clear unit_locations y_distance_matrix idx_row idx_col idx_good;

%% compute similarity
similarity_waveform = zeros(n_pairs, 1);
similarity_ISI = zeros(n_pairs, 1);
similarity_AutoCorr = zeros(n_pairs, 1);
similarity_PETH = zeros(n_pairs, 1);

% start parallel pool
if isempty(gcp('nocreate'))
    parpool();
end
disp('Start computing similarity!');
progBar = ProgressBar(n_pairs, ...
    'IsParallel', true, ...
    'Title', 'Computing Similarity', ...
    'UpdateRate', 1 ...
    );
progBar.setup([], [], []);

parfor k = 1:n_pairs
    idx_A = idx_unit_pairs(k,1);
    idx_B = idx_unit_pairs(k,2);
    
    if any(strcmpi(user_settings.motionEstimation.features, 'Waveform'))
        similarity_waveform(k) = waveformSimilarity(spikeInfo(idx_A), spikeInfo(idx_B),...
                user_settings.waveformCorrection.n_nearest_channels,...
                user_settings.motionEstimation.interpolate_algorithm);
    end
    
    if any(strcmpi(user_settings.motionEstimation.features, 'ISI'))
        similarity_ISI(k) = ISI_Similarity(spikeInfo(idx_A), spikeInfo(idx_B));
    end
    
    if any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr'))
        similarity_AutoCorr(k) = autocorrelogramSimilarity(spikeInfo(idx_A), spikeInfo(idx_B));
    end
    
    if any(strcmpi(user_settings.motionEstimation.features, 'PETH'))
        similarity_PETH(k) = PETH_Similarity(spikeInfo(idx_A), spikeInfo(idx_B));
    end
    
    updateParallel(1);
end
progBar.release();

fprintf('Computing similarity done! Saved to %s ...\n', fullfile(user_settings.output_folder, 'SimilarityForCorretion.mat'));
toc;

% save the similarity
if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'SimilarityForCorretion.mat'),...
        'similarity_waveform', 'similarity_ISI', 'similarity_AutoCorr', 'similarity_PETH', 'idx_unit_pairs', '-nocompression');
end

%% Pre-clustering
n_unit = length(spikeInfo);
sessions = [spikeInfo.SessionIndex];
n_session = max(sessions);
names_all = {'Waveform', 'ISI', 'AutoCorr', 'PETH'};
similarity_all = [similarity_waveform, similarity_ISI, similarity_AutoCorr, similarity_PETH];

similarity_names = user_settings.motionEstimation.features';
idx_names = zeros(1, length(similarity_names));
for k = 1:length(similarity_names)
    idx_names(k) = find(strcmpi(names_all, similarity_names{k}));
end
similarity_all = similarity_all(:, idx_names);

weights = ones(1, length(similarity_names))./length(similarity_names);
mean_similarity = sum(similarity_all.*weights, 2);

similarity_matrix = zeros(n_unit);
for k = 1:size(idx_unit_pairs, 1)
    similarity_matrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = mean_similarity(k);
    similarity_matrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = mean_similarity(k);
end
similarity_matrix(eye(size(similarity_matrix)) == 1) = 5;
%
for iter = 1:user_settings.motionEstimation.n_iter
    fprintf('Iteration %d starts!\n', iter);

    % HDBSCAN
    distance_matrix = 1./(1 + tanh(similarity_matrix));
    
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
    
    system([fullfile(user_settings.path_to_python), ' ',...
        fullfile(path_kilomatch, 'Functions/main_hdbscan.py'), ' ',...
        fullfile(user_settings.output_folder, 'HDBSCAN_settings.json')]);
    
    idx_cluster_hdbscan = double(readNPY(fullfile(user_settings.output_folder, 'ClusterIndices.npy')));
    % MATLAB starts from 1
    idx_cluster_hdbscan(idx_cluster_hdbscan >= 0) = idx_cluster_hdbscan(idx_cluster_hdbscan >= 0)+1;
    
    n_cluster = max(idx_cluster_hdbscan);
    hdbscan_matrix = zeros(size(similarity_matrix));
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
    is_matched = arrayfun(@(x)hdbscan_matrix(idx_unit_pairs(x,1), idx_unit_pairs(x,2)), 1:length(mean_similarity));

    if iter ~= user_settings.motionEstimation.n_iter        
        % LDA and update weights
        mdl = fitcdiscr(similarity_all, is_matched);
        temp = (mdl.Coeffs(1,2).Linear)';
        weights = temp./sum(temp);
        disp('Weights:');
        disp(strjoin(similarity_names, '   '));
        disp(weights);
        
        % update the similarity matrix
        mean_similarity = sum(similarity_all.*weights, 2);
        similarity_matrix = zeros(n_unit);
        for k = 1:size(idx_unit_pairs, 1)
            similarity_matrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = mean_similarity(k);
            similarity_matrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = mean_similarity(k);
        end
        similarity_matrix(eye(size(similarity_matrix)) == 1) = 5;
    end
end

% set the threshold based on LDA results
similarity_thres = mdl.Coeffs(1,2).Const ./ (-mdl.Coeffs(1,2).Linear(1)) .* weights(1);

similarity = sum(similarity_all.*weights, 2);
idx_good = find(similarity' > similarity_thres & is_matched == 1);
n_pairs_included = length(idx_good);

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
dx = zeros(1, length(idx_good));
idx_1 = zeros(1, length(idx_good));
idx_2 = zeros(1, length(idx_good));
for k = 1:length(idx_good)
    unit1 = idx_unit_pairs(idx_good(k), 1);
    unit2 = idx_unit_pairs(idx_good(k), 2);
    d_this = mean([spikeInfo(unit2).Location(2), spikeInfo(unit1).Location(2)]);

    idx_1(k) = session_pairs(idx_good(k),1);
    idx_2(k) = session_pairs(idx_good(k),2);
    dx(k) = spikeInfo(unit2).Location(2) - spikeInfo(unit1).Location(2);
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
    dx_block = dx(idx_block == k);
    idx_1_block = idx_1(idx_block == k);
    idx_2_block = idx_2(idx_block == k);

    if length(unique([idx_1_block, idx_2_block])) ~= n_session
        disp('Some sessions are not included! Motion estimation failed!');
        continue
    end

    options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8);
    loss_fun = @(x) sum((dx_block - (x(idx_2_block) - x(idx_1_block))).^2, 'all');

    p = fminunc(loss_fun, rand(1, n_session), options);
    p = p - mean(p);
    positions(k,:) = p;

    options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8, 'Display', 'none');
    p_boot = zeros(n_boot, n_session);

    progBar = ProgressBar(...
        n_boot, ...
        'Title', 'Computing 95CI' ...
        );
    for j = 1:n_boot
        idx_rand = randi(length(dx_block), 1, length(dx_block));
        dx_this = dx_block(idx_rand);
        idx_1_this = idx_1_block(idx_rand);
        idx_2_this = idx_2_block(idx_rand);
    
        loss_fun = @(x) sum((dx_this - (x(idx_2_this) - x(idx_1_this))).^2, 'all');
    
        p_this = fminunc(loss_fun, rand(1, n_session), options);
        p_boot(j,:) = p_this - mean(p_this);
    
        progBar([], [], []);
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
p_all_old = positions;
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
dx_left = dx;
for k = 1:nblock
    dx_left(idx_block==k) = dx(idx_block==k) - (positions(k, idx_2(idx_block==k)) - positions(k, idx_1(idx_block==k)));
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
            [min(dx), max(dx); min(dx), max(dx)]', 'shadedColor', colors_block(k,:));
    end
end

plot(ax_distance, depth, dx, '.', 'Color', colors(1,:));
plot(ax_distance, depth, dx_left, '.', 'Color', colors(2,:));

xlabel(ax_distance, 'Depth (um)');
ylabel(ax_distance, 'Residues (um)');

xlim(ax_distance, [depth_edges(1), depth_edges(end)]);
ylim(ax_distance, [min(dx), max(dx)]);

EasyPlot.cropFigure(fig);
if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Motion'));
end

% save data
if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'Motion.mat'), 'positions', 'depth_bins', 'nblock', '-nocompression');
end

% clear temp variables
clear similarity similarity_all similarity_matrix similarity_waveform similarity_ISI similarity_AutoCorr similarity_PETH;
