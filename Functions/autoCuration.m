function [hdbscan_matrix_curated, idx_cluster_hdbscan_curated, curation_pairs, curation_types, curation_type_names] = autoCuration(...
    user_settings, hdbscan_matrix, idx_cluster_hdbscan, good_matches_matrix, ...
    sessions, similarity_matrix, leafOrder)

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

fprintf('%d deleting steps are done!\n', num_before-num_after);
fprintf('%d clusters and %d pairs after removing bad units!\n',...
    n_cluster, (sum(hdbscan_matrix(:)) - n_unit)/2);

% merge when two or more adjacent clusters are similar and do not contain units from the same sessions
if user_settings.autoCuration.auto_merge
    % compute the location of each clusters
    hdbscan_matrix_raw = hdbscan_matrix;
    cluster_centers = zeros(1,n_cluster);
    for k = 1:n_cluster
        units = find(idx_cluster_hdbscan == k);
        temp = arrayfun(@(x)find(leafOrder == x), units);
        cluster_centers(k) = median(temp);
    end
    [~, cluster_id_sorted] = sort(cluster_centers);
    
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
    
            if ~any(good_matches_matrix(units, units_next) > 0, 'all')
                continue
            end
            
            similar_pairs = [similar_pairs; k, k+1, median(similarity_matrix(units, units_next), 'all')];
        end
        
        if ~isempty(similar_pairs)
            fprintf('Found %d possible merges!\n', size(similar_pairs, 1));
        end
    
        % merging
        for k = 1:size(similar_pairs, 1)
            if k < size(similar_pairs, 1) &&...
                    similar_pairs(k+1, 1) == similar_pairs(k,2) &&...
                    similar_pairs(k+1, 3) > similar_pairs(k,3)
                continue
            end
            
            flag = true;
            
            id = cluster_id_sorted(similar_pairs(k,1));
            id_next = cluster_id_sorted(similar_pairs(k,2));

            % update curation pairs
            units1 = find(idx_cluster_hdbscan == id);
            units2 = find(idx_cluster_hdbscan == id_next);
            for j = 1:length(units1)
                for i = 1:length(units2)
                    pair_this = sort([units1(j), units2(i)]);
                    curation_pairs = [curation_pairs; pair_this];
                    curation_types = [curation_types, 3];
                end
            end

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
    
    [num_same, num_before, num_after] = graphEditNumber(hdbscan_matrix_raw, hdbscan_matrix);
    assert(num_same == num_before);
    
    fprintf('%d merging steps are done!\n', -num_before+num_after);
    fprintf('%d clusters and %d pairs after merging good clusters!\n',...
        n_cluster, (sum(hdbscan_matrix(:)) - n_unit)/2);
    
    % find possible pairings for unpaired units
    disp('Checking the unpaired units!');
    
    count_merges = 0;
    
    n_cluster_new = n_cluster;
    flag_match = true;
    
    while flag_match
        flag_match =  false;
        idx_unpaired = find(idx_cluster_hdbscan == -1);
        for k = 1:length(idx_unpaired)
            unit = idx_unpaired(k);
            session_this = sessions(unit);
        
            idx_match = find(good_matches_matrix(unit,:) == 1);
    
            if isempty(idx_match)
                continue
            end
            
            [~, temp] = max(similarity_matrix(k, idx_match));
            idx_match = idx_match(temp);
            idx_cluster_new = idx_cluster_hdbscan(idx_match);
            
            if idx_cluster_new == -1
                continue
            end
            if any(session_this == sessions(idx_cluster_hdbscan == idx_cluster_new))
                continue
            end
            
            flag_match = true;
            count_merges = count_merges+1;

            % update curation pairs
            units1 = find(idx_cluster_hdbscan == idx_cluster_new);
            for j = 1:length(units1)
                pair_this = sort([units1(j), unit]);
                curation_pairs = [curation_pairs; pair_this];
                curation_types = [curation_types, 4];
            end

            idx_cluster_hdbscan(unit) = idx_cluster_new;
        end
    end
    n_cluster = n_cluster_new;
    fprintf('Merged %d unpaired units!\n', count_merges);
    
    % update hdbscan matrix
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
    
    fprintf('%d clusters and %d pairs after merging good unpaired units!\n',...
        n_cluster, (sum(hdbscan_matrix(:)) - n_unit)/2);
end    

hdbscan_matrix_curated = hdbscan_matrix;
idx_cluster_hdbscan_curated = idx_cluster_hdbscan;

% Plot the final results
fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, 1, 2,...
    'Width', 12,...
    'Height', 12,...
    'MarginBottom', 1,...
    'MarginLeft', 1,...
    'MarginRight', 0.2);

imagesc(ax_all{1}, hdbscan_matrix_curated(leafOrder,leafOrder));
imagesc(ax_all{2}, similarity_matrix(leafOrder,leafOrder));

EasyPlot.setCLim(ax_all{2}, [0, 4]);
EasyPlot.colorbar(ax_all{2},...
    'label', 'Similarity',...
    'MarginRight', 1);

EasyPlot.setXLim(ax_all, [0.5, length(leafOrder)+0.5]);
EasyPlot.setYLim(ax_all, [0.5, length(leafOrder)+0.5]);
title(ax_all{1}, 'Curated result');
title(ax_all{2}, 'Similarity matrix');

linkaxes([ax_all{1}, ax_all{2}]);

EasyPlot.cropFigure(fig);

if user_settings.save_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/CuratedResult'));
end


end