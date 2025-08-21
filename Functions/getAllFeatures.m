function [ISI_features, AutoCorr_features, PETH_features] = getAllFeatures(spikeInfo)
% GETALLFEATURES  Extract ISI, autocorrelogram, and PETH feature matrices.
%
% Gathers available features from a spikeInfo struct array:
% 1. Concatenates ISI distributions across units if .ISI exists.
% 2. Concatenates autocorrelograms across units if .AutoCorr exists.
% 3. Concatenates peristimulus time histograms across units if .PETH exists.
%
% Inputs:
%   spikeInfo           struct array (n_unit × 1)
%       Fields (optional):
%         .ISI            vector or matrix of ISI features
%         .AutoCorr       vector or matrix of autocorrelogram features
%         .PETH           vector or matrix of PETH features
%
% Outputs:
%   ISI_features        numeric matrix (n_unit × n_ISI_bins) or empty
%   AutoCorr_features   numeric matrix (n_unit × n_bins) or empty
%   PETH_features       numeric matrix (n_unit × n_timepoints) or empty
%
% Date:    20250821  
% Author:  Yue Huang 

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