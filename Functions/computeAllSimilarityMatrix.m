function [similarity_matrix_all, feature_names_all] = computeAllSimilarityMatrix( ...
    user_settings, waveforms_all, channel_locations, ISI_features, AutoCorr_features, PETH_features, feature_names, idx_unit_pairs)
% waveforms_all: n_unit x n_channel x n_sample array
% channel_locations: n_channel x 2 array
% ISI_features: n_unit x n_sample_ISI array
% AutoCorr_features: n_unit x n_sample_AutoCorr array
% PETH_features: n_unit x n_sample_PETH array
% feature_names: n_feature x 1 cell array

n_unit = size(waveforms_all, 1);

waveform_similarity_matrix = zeros(n_unit);
if any(strcmpi(feature_names, 'Waveform'))
    waveform_similarity_matrix = computeWaveformSimilarityMatrix(user_settings, waveforms_all, channel_locations);
end

% compute other similarity
ISI_similarity_matrix = zeros(n_unit);
if any(strcmpi(feature_names, 'ISI'))
    ISI_similarity_matrix = corrcoef(ISI_features');
    ISI_similarity_matrix(isnan(ISI_similarity_matrix)) = 0;
    ISI_similarity_matrix = atanh(ISI_similarity_matrix);
end

AutoCorr_similarity_matrix = zeros(n_unit);
if any(strcmpi(feature_names, 'AutoCorr'))
    AutoCorr_similarity_matrix = corrcoef(AutoCorr_features');
    AutoCorr_similarity_matrix(isnan(AutoCorr_similarity_matrix)) = 0;
    AutoCorr_similarity_matrix = atanh(AutoCorr_similarity_matrix);
end

PETH_similarity_matrix = zeros(n_unit);
if any(strcmpi(feature_names, 'PETH'))
    PETH_similarity_matrix = corrcoef(PETH_features');
    PETH_similarity_matrix(isnan(PETH_similarity_matrix)) = 0;
    PETH_similarity_matrix = atanh(PETH_similarity_matrix);
end

feature_names_all = {'Waveform', 'ISI', 'AutoCorr', 'PETH'};
similarity_matrix_all = cat(3,...
    waveform_similarity_matrix,...
    ISI_similarity_matrix,...
    AutoCorr_similarity_matrix,...
    PETH_similarity_matrix);

% plot the distribution of similarity
fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, 1, length(feature_names),...
    'Width', 5,...
    'Height', 5,...
    'MarginBottom', 1,...
    'MarginLeft', 1);

idx_pairs_in_matrix = sub2ind([n_unit, n_unit], idx_unit_pairs(:,1), idx_unit_pairs(:,2));

for k = 1:length(ax_all)
    idx_this = strcmpi(feature_names_all, feature_names{k});

    similarity_matrix_this = similarity_matrix_all(:, :, idx_this);
    histogram(ax_all{k}, similarity_matrix_this(idx_pairs_in_matrix), 'BinWidth', 0.2, 'Normalization', 'probability');
    xlabel(ax_all{k}, [feature_names{k}, ' similarity']);
    ylabel(ax_all{k}, 'Prob.');    
end

EasyPlot.setYLim(ax_all);
EasyPlot.cropFigure(fig);

if user_settings.save_intermediate_figures
    EasyPlot.exportFigure(fig, fullfile(user_settings.output_folder, 'Figures/AllSimilarity'));
end

end
