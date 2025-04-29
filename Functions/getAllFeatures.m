function [ISI_features, AutoCorr_features, PETH_features] = getAllFeatures(user_settings, spikeInfo)
% Get all features
if isfield(spikeInfo, 'ISI') &&...
        (any(strcmpi(user_settings.motionEstimation.features, 'ISI')) ||...
        any(strcmpi(user_settings.clustering.features, 'ISI')))
    ISI_features = cat(1, spikeInfo.ISI);
else
    ISI_features = [];
end

if isfield(spikeInfo, 'AutoCorr') &&...
        (any(strcmpi(user_settings.motionEstimation.features, 'AutoCorr')) ||...
        any(strcmpi(user_settings.clustering.features, 'AutoCorr')))
    AutoCorr_features = cat(1, spikeInfo.AutoCorr);
else
    AutoCorr_features = [];
end

if isfield(spikeInfo, 'PETH') &&...
        (any(strcmpi(user_settings.motionEstimation.features, 'PETH')) ||...
        any(strcmpi(user_settings.clustering.features, 'PETH')))
    PETH_features = cat(1, spikeInfo.PETH);
else
    PETH_features = [];
end

end