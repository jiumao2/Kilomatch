function Motion = computeMotion(user_settings, similarity_matrix, hdbscan_matrix, idx_unit_pairs, similarity_thres, sessions, locations)
idx_out = sub2ind(size(similarity_matrix), idx_unit_pairs(:,1), idx_unit_pairs(:,2));
good_matrix = (similarity_matrix > similarity_thres) & (hdbscan_matrix == 1);
idx_good = find(good_matrix(idx_out) == 1);

fprintf('%d pairs of units are included for drift estimation!\n', length(idx_good));

n_session = max(sessions);
session_pairs = sessions(idx_unit_pairs);

% get all the good pairs and their distance
session_pairs_good = session_pairs(idx_good,:)';

depth = zeros(1, length(idx_good));
dy = zeros(1, length(idx_good));
for k = 1:length(idx_good)
    unit1 = idx_unit_pairs(idx_good(k), 1);
    unit2 = idx_unit_pairs(idx_good(k), 2);
    depth(k) = mean(locations([unit1, unit2], 2));
    dy(k) = locations(unit2, 2) - locations(unit1, 2);
end

% compute the motion and 95CI
n_boot = 100;
linear_scale = 1/1000;

if length(unique(session_pairs_good)) ~= n_session
    disp('Some sessions are not included! Motion estimation failed!');
end

options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8);

Motion = struct('Linear', zeros(1, n_session), 'Constant', zeros(1, n_session), 'LinearScale', linear_scale);
if user_settings.waveformCorrection.linear_correction
    loss_fun = @(params)lossFunLinearCorrection(params, session_pairs_good, dy, depth, linear_scale, n_session);
    params = fminunc(loss_fun, rand(1, n_session*2-2), options);
    Motion.Linear = [0, params(1:n_session-1)];
    Motion.Constant = [0, params(n_session:end)];

    fprintf('Distortion range: %.2f%% ~ %.2f%% relative to the probe\n',...
        min(Motion.Linear)/range(depth)*100,...
        max(Motion.Linear)/range(depth)*100);

    mean_motion = Motion.LinearScale*Motion.Linear*mean(depth) + Motion.Constant;
    Motion.Constant = Motion.Constant - mean(mean_motion);

    mean_motion = Motion.LinearScale*Motion.Linear*mean(depth) + Motion.Constant;
    min_motion = Motion.LinearScale*Motion.Linear*min(depth) + Motion.Constant;
    max_motion = Motion.LinearScale*Motion.Linear*max(depth) + Motion.Constant;
else
    loss_fun = @(params)lossFunLinearDefault(params, session_pairs_good, dy);
    params = fminunc(loss_fun, rand(1, n_session), options);
    
    Motion.Constant = params - mean(params);

    mean_motion = Motion.Constant;
    min_motion = [];
    max_motion = [];
end

options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8, 'Display', 'none');
mean_motion_boot = zeros(n_boot, n_session);
min_motion_boot = zeros(n_boot, n_session);
max_motion_boot = zeros(n_boot, n_session);

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
    idx_rand = randi(length(dy), 1, length(dy));
    dy_this = dy(idx_rand);
    session_pairs_this = session_pairs_good(:, idx_rand);
    depth_this = depth(idx_rand);

    if user_settings.waveformCorrection.linear_correction
        loss_fun = @(params)lossFunLinearCorrection(params, session_pairs_good, dy, depth, linear_scale, n_session);
        params = fminunc(loss_fun, rand(1, n_session*2-2), options);
        p_linear = [0, params(1:n_session-1)];
        p_constant = [0, params(n_session:end)];
    
        mean_motion_boot(j,:) = linear_scale*p_linear*mean(depth) + p_constant;
        min_motion_boot(j,:) = linear_scale*p_linear*min(depth) + p_constant;
        max_motion_boot(j,:) = linear_scale*p_linear*max(depth) + p_constant;
    else
        loss_fun = @(params)lossFunLinearDefault(params, session_pairs_good, dy);
        params = fminunc(loss_fun, rand(1, n_session), options);
        mean_motion_boot(j,:) = params - mean(params);
    end

    updateParallel(1);
end
progBar.release();

mean_motion_ci95 = zeros(2, n_session);
min_motion_ci95 = zeros(2, n_session);
max_motion_ci95 = zeros(2, n_session);
for j = 1:n_session
    mean_motion_ci95(1,j) = prctile(mean_motion_boot(:,j), 2.5);
    mean_motion_ci95(2,j) = prctile(mean_motion_boot(:,j), 97.5);

    min_motion_ci95(1,j) = prctile(min_motion_boot(:,j), 2.5);
    min_motion_ci95(2,j) = prctile(min_motion_boot(:,j), 97.5);

    max_motion_ci95(1,j) = prctile(max_motion_boot(:,j), 2.5);
    max_motion_ci95(2,j) = prctile(max_motion_boot(:,j), 97.5);
end    

% plot the final similarity score distribution
similarity = similarity_matrix(idx_out);

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
title(ax_similarity, [num2str(length(idx_good)), ' pairs are included']);
EasyPlot.cropFigure(fig);

if user_settings.save_intermediate_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/SimilarityThresholdForCorrection'), 'dpi', 300);
end

% Plot the drift
% compute the residues
location_estimated = depth*Motion.LinearScale.*Motion.Linear(session_pairs_good) + Motion.Constant(session_pairs_good);
dy_estimated = location_estimated(2,:) - location_estimated(1,:);

dy_left = dy - dy_estimated;

%% Make the plot
fig = EasyPlot.figure();
ax_drift = EasyPlot.axes(fig, 'Height', 3, 'Width', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 4);

ax_distance = EasyPlot.createAxesAgainstAxes(fig, ax_drift, 'right',...
    'Width', 8);

colors = [0.7,0.7,0.7; lines(1)];

% drift
h1 = plot(ax_drift, 1:n_session, mean_motion, '-', 'Color', 'k');
EasyPlot.plotShaded(ax_drift, 1:n_session, mean_motion_ci95, 'shadedColor', 'k');

if user_settings.waveformCorrection.linear_correction
    h2 = plot(ax_drift, 1:n_session, min_motion, '-', 'Color', 'b');
    EasyPlot.plotShaded(ax_drift, 1:n_session, min_motion_ci95, 'shadedColor', 'b');
    h3 = plot(ax_drift, 1:n_session, max_motion, '-', 'Color', 'r');
    EasyPlot.plotShaded(ax_drift, 1:n_session, max_motion_ci95, 'shadedColor', 'r');

    h_legend = EasyPlot.legend(ax_drift, {'Mean motion', 'Bottom motion', 'Top motion'},...
        'location', 'northeastoutside',...
        'selectedPlots', [h1,h2,h3],...
        'lineLength', 0.5);
    EasyPlot.move(h_legend, 'dx', -1);
end

xlabel(ax_drift, 'Session');
ylabel(ax_drift, 'Motion (μm)');

% raw distance
plot(ax_distance, depth, dy, '.', 'Color', colors(1,:));
plot(ax_distance, depth, dy_left, '.', 'Color', colors(2,:));

xlabel(ax_distance, 'Depth (μm)');
ylabel(ax_distance, 'Residues (μm)');

xlim(ax_distance, [min(depth), max(depth)]);
ylim(ax_distance, [min(dy), max(dy)]);

EasyPlot.cropFigure(fig);
if user_settings.save_intermediate_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Motion'));
end

% save data
if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'Motion.mat'), 'Motion', '-nocompression');
end

end