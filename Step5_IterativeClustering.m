%% compute similarity matrix
n_unit = length(spikeInfo);
sessions = [spikeInfo.SessionIndex];
n_session = max(sessions);
similarity_matrix = zeros(n_unit);

names_all = {'Waveform', 'ISI', 'AutoCorr', 'PETH'};
similarity_all = [similarity_waveform, similarity_ISI, similarity_AutoCorr, similarity_PETH];

similarity_names = user_settings.motionEstimation.features;
idx_names = zeros(1, length(similarity_names));
for k = 1:length(similarity_names)
    idx_names{k} = find(strcmpi(names_all, similarity_names{k}));
end
similarity_all = similarity_all(:, idx_names);

weights = ones(1, length(similarity_names))./length(similarity_names);
mean_similarity = sum(similarity_all.*weights, 2);

for k = 1:size(idx_unit_pairs, 1)
    similarity_matrix(idx_unit_pairs(k,1), idx_unit_pairs(k,2)) = mean_similarity(k);
    similarity_matrix(idx_unit_pairs(k,2), idx_unit_pairs(k,1)) = mean_similarity(k);
end
similarity_matrix(eye(size(similarity_matrix)) == 1) = 5;

for iter = 1:user_settings.clustering.n_iter
    fprintf('Iteration %d starts!\n', iter);

    % HDBSCAN
    distance_matrix = 1./(1 + tanh(similarity_matrix));
    
    HDBSCAB_settings.min_samples = 2; % The number of samples in a neighborhood for a point to be considered as a core point.
    % This includes the point itself. When None, defaults to min_cluster_size.
    HDBSCAB_settings.cluster_selection_epsilon = 0; % A distance threshold. 
    % Clusters below this value will be merged. This is the minimum epsilon allowed.
    
    HDBSCAB_settings.min_cluster_size = 2;
    HDBSCAB_settings.max_cluster_size = n_session;
    HDBSCAB_settings.metric = 'precomputed';
    HDBSCAB_settings.data_folder = user_settings.output_folder;
    
    json_text = jsonencode(HDBSCAB_settings);
    
    writelines(json_text, fullfile(user_settings.output_folder, 'HDBSCAN_settings.json'));
    writeNPY(distance_matrix, fullfile(user_settings.output_folder, 'DistanceMatrix.npy'));
    
    system([fullfile(user_settings.path_to_python), ' ',...
        fullfile(path_kilomatch, 'Functions/hdbscan.py'), ' ',...
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

    if iter ~= user_settings.clustering.n_iter
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

% do hierachical clustering on dissimilarity matrix
dissimilarity_matrix = 1-tanh(similarity_matrix);
dissimilarity_matrix(eye(size(dissimilarity_matrix)) == 1) = 0;
dissimilarity_vector = squareform(dissimilarity_matrix);
Z = linkage(dissimilarity_vector, 'average');
coph_corr = cophenet(Z, dissimilarity_vector);
disp(['Cophenetic Correlation Coefficient = ', num2str(coph_corr)]);

leafOrder = optimalleaforder(Z, dissimilarity_matrix, 'Transformation', @(x)1-x);

similarity = sum(similarity_all.*weights, 2);
edges = 0:0.01:5;
hist_matched = histcounts(similarity(is_matched == 1), edges, 'Normalization', 'probability');
hist_unmatched = histcounts(similarity(is_matched == 0), edges, 'Normalization', 'probability');
hist_diff = hist_matched - hist_unmatched;
idx_crossed = find(sign(hist_diff) > 0, 1);
thres = edges(idx_crossed);

good_matches_matrix = zeros(size(similarity_matrix));
idx_good_matches = find((sum(similarity_all.*weights, 2) - thres) > 0);
for k = 1:length(idx_good_matches)
    good_matches_matrix(idx_unit_pairs(idx_good_matches(k), 1), idx_unit_pairs(idx_good_matches(k), 2)) = 1;
    good_matches_matrix(idx_unit_pairs(idx_good_matches(k), 2), idx_unit_pairs(idx_good_matches(k), 1)) = 1;
end
good_matches_matrix(eye(size(good_matches_matrix)) == 1) = 1;
%% Save the results
save(fullfile(user_settings.output_folder, 'ClusteringResults.mat'),...
    'sessions', 'n_session', 'weights', 'similarity_names', 'similarity_all',...
    'thres', 'good_matches_matrix',...
    'similarity_matrix', 'distance_matrix', 'dissimilarity_matrix', 'leafOrder',...
    'HDBSCAB_settings', 'idx_cluster_hdbscan', 'hdbscan_matrix', 'n_cluster');

%% Plot the similarity histogram
fig = EasyPlot.figure();
ax = EasyPlot.axes(fig,...
    'Width', 5,...
    'Height', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

% weights = model_coeff_linear./sum(model_coeff_linear);
% thres = -model_coeff_const./sum(model_coeff_linear);
sim_this = sum(similarity_all.*weights, 2);
histogram(ax, sim_this(is_matched==0), 'Normalization', 'probability', 'BinWidth', 0.05, 'FaceColor', 'k');
histogram(ax, sim_this(is_matched==1), 'Normalization', 'probability', 'BinWidth', 0.05, 'FaceColor', 'b');
xline(ax, thres, 'k:', 'LineWidth', 2);
xlabel(ax, 'Similarity');
ylabel(ax, 'Probability');

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/SummedSimilarityDistribution'));
%% Plot the results
fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, 1, 3,...
    'Width', 12,...
    'Height', 12,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 0.2);

imagesc(ax_all{1}, hdbscan_matrix(leafOrder,leafOrder));
imagesc(ax_all{2}, good_matches_matrix(leafOrder,leafOrder));
imagesc(ax_all{3}, similarity_matrix(leafOrder,leafOrder));

EasyPlot.setCLim(ax_all{3}, [0, 4]);
h = EasyPlot.colorbar(ax_all{3},...
    'label', 'Similarity',...
    'MarginRight', 1);

EasyPlot.setXLim(ax_all, [0.5, length(leafOrder)+0.5]);
EasyPlot.setYLim(ax_all, [0.5, length(leafOrder)+0.5]);
title(ax_all{1}, 'IHDBSCAN');
title(ax_all{2}, 'Good matches');
title(ax_all{3}, 'Similarity matrix');

linkaxes([ax_all{1}, ax_all{2}, ax_all{3}]);

EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/ClusteringResult'));
savefig(fig, fullfile(user_settings.output_folder, 'Figures/ClusteringResult.fig'));


