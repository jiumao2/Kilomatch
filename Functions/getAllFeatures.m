function [ISI_features, AutoCorr_features, PETH_features] = getAllFeatures(user_settings, spikeInfo)
% Get all features
if isfield(spikeInfo, 'ISI')
    ISI_features = cat(1, spikeInfo.ISI);
else
    ISI_features = [];
end

if isfield(spikeInfo, 'AutoCorr')
    AutoCorr_features = cat(1, spikeInfo.AutoCorr);
else
    AutoCorr_features = [];
end

if isfield(spikeInfo, 'PETH')
    PETH_features = cat(1, spikeInfo.PETH);
else
    PETH_features = [];
end

end