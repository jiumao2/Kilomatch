function [hdbscan_matrix, idx_cluster_hdbscan]= autoCuration(hdbscan_matrix, idx_cluster_hdbscan, leafOrder,...
    similarity_matrix, sessions, user_settings)

merge_thres = user_settings.autoCuration.merge_threshold;
reject_thres = user_settings.autoCuration.reject_threshold;

% remove bad units in the cluster if any similarity == 0
n_cluster = max(idx_cluster_hdbscan);
fprintf('%d clusters and %d pairs before removing bad units!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - length(leafOrder))/2);

for k = 1:n_cluster
    units = find(idx_cluster_hdbscan == k);
    sessions_this = sessions(units);
    similarity_matrix_this = similarity_matrix(units, units);

    while length(sessions_this) ~= length(unique(sessions_this))...
            || any(similarity_matrix_this<reject_thres, 'all')
        idx_remove = [];
        for j = 1:length(sessions_this)
            for i = j+1:length(sessions_this)
                if sessions_this(i) == sessions_this(j) || similarity_matrix_this(i,j)<reject_thres
                    similarity_i = mean(similarity_matrix_this(i,:), 'all');
                    similarity_j = mean(similarity_matrix_this(j,:), 'all');

                    if similarity_i <= similarity_j
                        idx_remove = [idx_remove, i];
                    else
                        idx_remove = [idx_remove, j];
                    end
                end
            end
        end

        idx_cluster_hdbscan(units(idx_remove)) = -1;
        units(idx_remove) = [];
        sessions_this(idx_remove) = [];
        similarity_matrix_this(idx_remove, :) = [];
        similarity_matrix_this(:, idx_remove) = [];
    end
end

% update the clusters and hdbscan matrix
idx_remove = [];
for k = 1:n_cluster
    units = find(idx_cluster_hdbscan == k);
    if length(units) <= 1
        % remove the cluster
        idx_cluster_hdbscan(units) = -1;
        idx_remove = [idx_remove, k];
    end
end

for k = length(idx_remove):-1:1
    idx_this = idx_remove(k);
    idx_cluster_hdbscan(idx_cluster_hdbscan>=idx_this) = idx_cluster_hdbscan(idx_cluster_hdbscan>=idx_this)-1;
end

assert(length(unique(idx_cluster_hdbscan)) == max(idx_cluster_hdbscan)+1);

% update matrix
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
hdbscan_matrix(eye(length(leafOrder)) == 1) = 1;

fprintf('%d clusters and %d pairs after removing bad units!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - length(leafOrder))/2);

% merge when two or more adjacent clusters are similar and do not contain units from the same sessions

% compute the location of each clusters
cluster_centers = zeros(1,n_cluster);
for k = 1:n_cluster
    units = find(idx_cluster_hdbscan == k);
    temp = arrayfun(@(x)find(leafOrder == x), units);
    cluster_centers(k) = median(temp);
end
[~, cluster_id_sorted] = sort(cluster_centers);

similarity_between_clusters = [];
flag = true;
while flag
    similar_pairs = []; % idA, idB, similarity
    flag = false;
    for k = 1:length(cluster_id_sorted)-1
        id = cluster_id_sorted(k);
        units = find(idx_cluster_hdbscan == id);
        sessions_this = sessions(units);

        id_next = cluster_id_sorted(k+1);
        units_next = find(idx_cluster_hdbscan == id_next);
        sessions_next = sessions(units_next);
        if ~isempty(intersect(sessions_this, sessions_next))
            continue
        end
        
        if any(similarity_matrix(units, units_next) < reject_thres, 'all')
            continue
        end
        
        similar_pairs = [similar_pairs; k, k+1, median(similarity_matrix(units, units_next), 'all')];
    end
    
    if ~isempty(similar_pairs)
        fprintf('Found %d possible merges!\n', sum(similar_pairs(:,3) > merge_thres));
    end

    if isempty(similarity_between_clusters)
        similarity_between_clusters = similar_pairs(:,3);
    end

    % merging
    for k = 1:size(similar_pairs, 1)
        if similar_pairs(k,3) < merge_thres
            continue
        end

        if k < size(similar_pairs, 1) &&...
                similar_pairs(k+1, 1) == similar_pairs(k,2) &&...
                similar_pairs(k+1, 3) > similar_pairs(k,3)
            continue
        end
        
        flag = true;

        id = cluster_id_sorted(similar_pairs(k,1));
        id_next = cluster_id_sorted(similar_pairs(k,2));
        idx_cluster_hdbscan(idx_cluster_hdbscan == id_next) = id;
    end

    % update cluster info
    max_id = max(idx_cluster_hdbscan);
    for k = max_id:-1:1
        units = find(idx_cluster_hdbscan == k);
        if isempty(units)
            idx_cluster_hdbscan(idx_cluster_hdbscan>=k) = idx_cluster_hdbscan(idx_cluster_hdbscan>=k)-1;
        end
    end

    n_cluster = max(idx_cluster_hdbscan);
    assert(length(unique(idx_cluster_hdbscan)) == n_cluster+1);

    % update cluster centers
    cluster_centers = zeros(1,n_cluster);
    for k = 1:n_cluster
        units = find(idx_cluster_hdbscan == k);
        temp = arrayfun(@(x)find(leafOrder == x), units);
        cluster_centers(k) = median(temp);
    end
    [~, cluster_id_sorted] = sort(cluster_centers);
end

% update hdbscan matrix
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
hdbscan_matrix(eye(length(leafOrder)) == 1) = 1;

fprintf('%d clusters and %d pairs after merging good clusters!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - length(leafOrder))/2);

%% find possible pairings for unpaired units
disp('Checking the unpaired units!');
count_merges = 0;
similarity_unpaired_units = [];
for k = 1:size(similarity_matrix, 1)
    if idx_cluster_hdbscan(k) ~= -1
        continue
    end

    center_this = find(leafOrder == k);
    session_this = sessions(k);

    % find the left and right nearest cluster
    idx_right = find(cluster_centers - center_this > 0);
    idx_left = find(cluster_centers - center_this <= 0);
    [~, idxA] = min(abs(center_this - cluster_centers(idx_left)));
    [~, idxB] = min(abs(center_this - cluster_centers(idx_right)));

    cluster_ids = [idx_left(idxA), idx_right(idxB)];
    similarity_this = zeros(1, length(cluster_ids));
    for j = 1:length(cluster_ids)
        units = find(idx_cluster_hdbscan == cluster_ids(j));

        if any(sessions(units) == session_this)
            continue
        end

        if any(similarity_matrix(k, units) < reject_thres)
            continue
        end

        similarity_this(j) = mean(similarity_matrix(k, units));
        similarity_unpaired_units = [similarity_unpaired_units, similarity_this(j)];
    end

    [~, idx_max] = max(similarity_this);
    if similarity_this(idx_max) > merge_thres
        idx_cluster_hdbscan(k) = cluster_ids(idx_max);
        count_merges = count_merges+1;
    end
end

fprintf('Merged %d unpaired units!\n', count_merges);

fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, 1, 2,...
    'Width', 6,...
    'Height', 6,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

histogram(ax_all{1}, similarity_between_clusters, 'BinWidth', 0.1);
xline(ax_all{1}, merge_thres, 'k:', 'LineWidth', 2);
xlabel(ax_all{1}, 'Similarity between clusters');
ylabel(ax_all{1}, 'Number of pairs');

histogram(ax_all{2}, similarity_unpaired_units, 'BinWidth', 0.1);
xline(ax_all{2}, merge_thres, 'k:', 'LineWidth', 2);
xlabel(ax_all{2}, 'Similarity between unpaired units and clusters');
ylabel(ax_all{2}, 'Number of pairs');
EasyPlot.cropFigure(fig);
EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/CurationSimilarityDistribution'));

close all;

% update hdbscan matrix
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
hdbscan_matrix(eye(length(leafOrder)) == 1) = 1;

fprintf('%d clusters and %d pairs after merging good unpaired units!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - length(leafOrder))/2);
end



