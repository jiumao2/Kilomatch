function waveform_similarity_matrix = computeWaveformSimilarityMatrix(user_settings, waveforms_all, channel_locations)
% waveforms_all: n_unit x n_channel x n_sample array
% channel_locations: n_channel x 2 array
% sessions: n_unit x 1 array
% locations: n_unit x 2 array

n_nearest_channels = user_settings.waveformCorrection.n_nearest_channels;
n_unit = size(waveforms_all, 1);
n_templates = size(waveforms_all, 4);

idx_nearest = knnsearch(channel_locations, channel_locations, 'K', n_nearest_channels);
idx_nearest_sorted = sort(idx_nearest, 2);

[idx_nearest_unique, ~, idx_groups] = unique(idx_nearest_sorted, 'rows');

waveform_similarity_matrix = zeros(n_unit, n_unit, n_templates);
for i_template = 1:n_templates
    waveform_similarity_matrix_this = zeros(n_unit);
    ptt = max(waveforms_all(:,:,:,i_template),[],3) - min(waveforms_all(:,:,:,i_template), [], 3);
    [~, ch] = max(ptt, [], 2);
    
    progBar = ProgressBar(size(idx_nearest_unique, 1), ...
        'Title', 'Computing waveform features', ...
        'UpdateRate', 1);
    for k = 1:size(idx_nearest_unique, 1)
        idx_included = find(idx_groups == k);
        idx_units = find(arrayfun(@(x)any(idx_included==x), ch));
    
        if isempty(idx_units)
            progBar([], [], []);
            continue
        end
    
        waveform_this = reshape(waveforms_all(:, idx_nearest_unique(k,:), :, i_template), n_unit, [])';

        temp = corr(waveform_this(:, idx_units), waveform_this);
        temp(isnan(temp)) = 0;
        temp = atanh(temp);
        
        waveform_similarity_matrix_this(idx_units,:) = temp;
    
        progBar([], [], []);
    end
    progBar.release();
    
    waveform_similarity_matrix_this = max(...
        cat(3, waveform_similarity_matrix_this, waveform_similarity_matrix_this'),...
        [], 3);

    waveform_similarity_matrix(:,:,i_template) = waveform_similarity_matrix_this;
end

waveform_similarity_matrix = max(waveform_similarity_matrix, [], 3);

end