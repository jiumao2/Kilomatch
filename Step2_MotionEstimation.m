% Estimating the motion of the electrode
disp('---------------Motion Estimation---------------');

% Reserve memory for the variables
similarity = zeros(1e7, 1);
idx_unit_pairs = zeros(1e7, 2);

count = 0;
% find the best unit pairs between days and the corresponding movement
for k = 1:length(spikeInfo)
    for j = k+1:length(spikeInfo)
        % excluded unit pairs from the same session
        if spikeInfo(k).SessionIndex ==spikeInfo(j).SessionIndex
            continue
        end
        
        % exclude unit pairs that are far away
        if abs(spikeInfo(k).Location(2) - spikeInfo(j).Location(2)) > user_settings.motionEstimation.max_motion_distance
            continue
        end
        
        count = count+1;

        similarity(count) = waveformSimilarity(spikeInfo(k), spikeInfo(j),...
            user_settings.waveformCorrection.n_nearest_channels,...
            user_settings.waveformCorrection.interpolate_algorithm);

        idx_unit_pairs(count,:) = [k,j];
    end

    if mod(k, 20) == 1
        toc;
        fprintf('%d / %d done!\n', k, length(spikeInfo));
    end
end

similarity = similarity(1:count);
idx_unit_pairs = idx_unit_pairs(1:count,:);

% save the similarity
save(fullfile(user_settings.output_folder, 'SimilarityForCorretion.mat'), 'similarity', 'idx_unit_pairs');

similarity_thres = user_settings.motionEstimation.similarity_threshold;
nblock = user_settings.motionEstimation.n_block;

n_pairs_included = sum(similarity>similarity_thres);

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
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/SimilarityThresholdForCorrection'));

fprintf('%d pairs of units are included for drift estimation!\n', sum(similarity>similarity_thres));

%% compute drift
n_session = max([spikeInfo.SessionIndex]);
session_pairs = [[spikeInfo(idx_unit_pairs(:,1)).SessionIndex]', [spikeInfo(idx_unit_pairs(:,2)).SessionIndex]'];

% get all the good pairs and their distance
depth = [];
dx = [];
idx_1 = [];
idx_2 = [];
idx_good = find(similarity > similarity_thres);
for k = 1:length(idx_good)
    unit1 = idx_unit_pairs(idx_good(k), 1);
    unit2 = idx_unit_pairs(idx_good(k), 2);
    d_this = mean([spikeInfo(unit2).Location(2), spikeInfo(unit1).Location(2)]);

    idx_1 = [idx_1, session_pairs(idx_good(k),1)];
    idx_2 = [idx_2, session_pairs(idx_good(k),2)];
    dx = [dx, spikeInfo(unit2).Location(2) - spikeInfo(unit1).Location(2)];
    depth = [depth, d_this];
end

depth_edges = linspace(min(depth), max(depth), nblock+1);
depth_bins = 0.5*(depth_edges(1:end-1) + depth_edges(2:end));
idx_block = findNearestPoint(depth_bins, depth);

% compute the motion and 95CI
n_boot = 1000;
positions = NaN(nblock, n_session);
positions_ci95 = zeros(2, nblock, n_session);
for k = 1:nblock
    dx_block = dx(idx_block == k);
    idx_1_block = idx_1(idx_block == k);
    idx_2_block = idx_2(idx_block == k);

    if length(unique([idx_1_block, idx_2_block])) ~= n_session
        continue
    end

    options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8);
    loss_fun = @(x) sum((dx_block - (x(idx_2_block) - x(idx_1_block))).^2, 'all');

    p = fminunc(loss_fun, rand(1, n_session), options);
    p = p - mean(p);
    positions(k,:) = p;

    options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8, 'Display', 'none');
    p_boot = zeros(n_boot, n_session);
    for j = 1:n_boot
        idx_rand = randi(length(dx_block), 1, length(dx_block));
        dx_this = dx_block(idx_rand);
        idx_1_this = idx_1_block(idx_rand);
        idx_2_this = idx_2_block(idx_rand);
    
        loss_fun = @(x) sum((dx_this - (x(idx_2_this) - x(idx_1_this))).^2, 'all');
    
        p_this = fminunc(loss_fun, rand(1, n_session), options);
        p_boot(j,:) = p_this - mean(p_this);
    
        if mod(j, 100) == 1
            fprintf('%d / %d done!\n', j, n_boot);
        end
    end

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

% xline(ax_distance, d_edges, 'k:', 'LineWidth', 2);

xlabel(ax_distance, 'Depth (um)');
ylabel(ax_distance, 'Residues (um)');

xlim(ax_distance, [depth_edges(1), depth_edges(end)]);
ylim(ax_distance, [min(dx), max(dx)]);

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Motion'));

% save data
save(fullfile(user_settings.output_folder, 'Motion.mat'), 'positions', 'depth_bins', 'nblock');