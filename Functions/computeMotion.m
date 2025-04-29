function motion = computeMotion(user_settings, similarity_matrix, hdbscan_matrix, idx_unit_pairs, similarity_thres, sessions, locations)
idx_out = sub2ind(size(similarity_matrix), idx_unit_pairs(:,1), idx_unit_pairs(:,2));
good_matrix = (similarity_matrix > similarity_thres) & (hdbscan_matrix == 1);
idx_good = find(good_matrix(idx_out) == 1);

fprintf('%d pairs of units are included for drift estimation!\n', length(idx_good));

n_session = max(sessions);
session_pairs = sessions(idx_unit_pairs);

% get all the good pairs and their distance
depth = zeros(1, length(idx_good));
dy = zeros(1, length(idx_good));
idx_1 = zeros(1, length(idx_good));
idx_2 = zeros(1, length(idx_good));
for k = 1:length(idx_good)
    unit1 = idx_unit_pairs(idx_good(k), 1);
    unit2 = idx_unit_pairs(idx_good(k), 2);
    d_this = mean([locations(unit2, 2), locations(unit1, 2)]);

    idx_1(k) = session_pairs(idx_good(k),1);
    idx_2(k) = session_pairs(idx_good(k),2);
    dy(k) = locations(unit2, 2) - locations(unit1, 2);
    depth(k) = d_this;
end

% compute the motion and 95CI
n_boot = 100;

if length(unique([idx_1, idx_2])) ~= n_session
    disp('Some sessions are not included! Motion estimation failed!');
end

options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8);
loss_fun = @(y) sum((dy - (y(idx_2) - y(idx_1))).^2, 'all');

p = fminunc(loss_fun, rand(1, n_session), options);
motion = p - mean(p);

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
    idx_rand = randi(length(dy), 1, length(dy));
    dy_this = dy(idx_rand);
    idx_1_this = idx_1(idx_rand);
    idx_2_this = idx_2(idx_rand);

    loss_fun = @(y) sum((dy_this - (y(idx_2_this) - y(idx_1_this))).^2, 'all');

    p_this = fminunc(loss_fun, rand(1, n_session), options);
    p_boot(j,:) = p_this - mean(p_this);

    updateParallel(1);
end
progBar.release();

motion_ci95 = zeros(2, n_session);
for j = 1:n_session
    motion_ci95(1,j) = prctile(p_boot(:,j), 2.5);
    motion_ci95(2,j) = prctile(p_boot(:,j), 97.5);
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

if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/SimilarityThresholdForCorrection'));
end

% Plot the drift
% compute the residues
dy_left = dy - (motion(idx_2) - motion(idx_1));

fig = EasyPlot.figure();
ax_drift = EasyPlot.axes(fig, 'Height', 3, 'Width', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 1);

ax_distance = EasyPlot.createAxesAgainstAxes(fig, ax_drift, 'right',...
    'Width', 8);

colors = [0.7,0.7,0.7; lines(1)];

% drift
plot(ax_drift, 1:n_session, motion, '-', 'Color', 'k');
EasyPlot.plotShaded(ax_drift, 1:n_session, motion_ci95, 'shadedColor', 'k');

xlabel(ax_drift, 'Session');
ylabel(ax_drift, 'Motion (um)');

% raw distance
plot(ax_distance, depth, dy, '.', 'Color', colors(1,:));
plot(ax_distance, depth, dy_left, '.', 'Color', colors(2,:));

xlabel(ax_distance, 'Depth (um)');
ylabel(ax_distance, 'Residues (um)');

xlim(ax_distance, [min(depth), max(depth)]);
ylim(ax_distance, [min(dy), max(dy)]);

EasyPlot.cropFigure(fig);
if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/Motion'));
end

% save data
if user_settings.save_intermediate_results
    save(fullfile(user_settings.output_folder, 'Motion.mat'), 'motion', '-nocompression');
end

end