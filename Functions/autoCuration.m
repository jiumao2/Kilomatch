function [hdbscan_matrix, idx_cluster_hdbscan,...
    curation_pairs, curation_types, curation_type_names, num_removal] = autoCuration(...
    user_settings, hdbscan_matrix, idx_cluster_hdbscan, good_matches_matrix, ...
    sessions, similarity_matrix)

curation_type_names = {'Removal_SameSession', 'Removal_LowSimilarity', 'Merge_Cluster', 'Merge_Unit'};
curation_pairs = [];
curation_types = [];

n_unit = size(similarity_matrix, 1);

n_cluster = max(idx_cluster_hdbscan);
fprintf('%d clusters and %d pairs before removing bad units!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - size(hdbscan_matrix, 1))/2);

hdbscan_matrix_raw = hdbscan_matrix;

% remove bad units in the cluster if two or more units are originated from the same sessions 
% or any similarity < reject_thres
for k = 1:n_cluster
    units = find(idx_cluster_hdbscan == k);
    sessions_this = sessions(units);
    similarity_matrix_this = similarity_matrix(units, units);

    while length(sessions_this) ~= length(unique(sessions_this))
        idx_remove = [];
        for j = 1:length(sessions_this)
            for i = j+1:length(sessions_this)
                if sessions_this(i) == sessions_this(j)
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

        % update curation pairs
        for j = 1:length(idx_remove)
            unit1 = units(idx_remove(j));
            for i = 1:length(units)
                if any(idx_remove(1:j) == i)
                    continue
                end
                unit2 = units(i);
                pair_this = sort([unit1, unit2]);
                curation_pairs = [curation_pairs; pair_this];
                curation_types = [curation_types, 1];
            end
        end

        idx_cluster_hdbscan(units(idx_remove)) = -1;
        units(idx_remove) = [];
        sessions_this(idx_remove) = [];
        similarity_matrix_this(idx_remove, :) = [];
        similarity_matrix_this(:, idx_remove) = [];
    end
end

% split a cluster if there is a clear boundary in good_matches_matrix
if user_settings.autoCuration.auto_split
    n_cluster_new = n_cluster;
    for k = 1:n_cluster
        units = find(idx_cluster_hdbscan == k);
        graph_this = graph(good_matches_matrix(units, units));
        idx_sub_clusters = conncomp(graph_this);
    
        n_sub_clusters = max(idx_sub_clusters);
        if n_sub_clusters <= 1
            continue
        end

        for j = 2:n_sub_clusters
            units_this = units(idx_sub_clusters == j);
            idx_cluster_hdbscan(units_this) = n_cluster_new+j-1;

            % update curation pairs
            for i = 1:length(units_this)
                unit1 = units_this(i);
                for ii = 1:length(units)
                    if any(units_this == units(ii))
                        continue
                    end
                    unit2 = units(ii);
                    pair_this = sort([unit1, unit2]);
                    curation_pairs = [curation_pairs; pair_this];
                    curation_types = [curation_types, 2];
                end
            end
        end
    
        n_cluster_new = n_cluster_new + n_sub_clusters - 1;
    end
    n_cluster = n_cluster_new;
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
hdbscan_matrix(eye(n_unit) == 1) = 1;

[num_same, num_before, num_after] = graphEditNumber(hdbscan_matrix_raw, hdbscan_matrix);
assert(num_same == num_after);

num_removal = num_before-num_after;
fprintf('%d deleting steps are done!\n', num_removal);
fprintf('%d clusters and %d pairs after removing bad units!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - n_unit)/2);

end