function [auto_corr, lag] = computeAutoCorr(spike_times, window, binwidth)
% COMPUTEAUTOCORR  Compute the autocorrelogram using spike times without normalization.
%
% computeAutoCorr calculates unnormalized autocorrelation counts by computing
% time differences between all spike pairs, binning them within ±window
% milliseconds at a resolution of binwidth. The lag vector spans from
% –window to +window in steps of binwidth.
%
% Inputs:
%   spike_times      double (1 × n)  
%       Spike times in milliseconds.
%
%   window           double (1 × 1)  
%       Half‐width of the correlogram window in ms (default: 300).
%
%   binwidth         double (1 × 1)  
%       Bin width for time differences in ms (default: 1).
%
% Outputs:
%   auto_corr        double (1 × (2*window+1))  
%       Autocorrelation counts for each lag bin.
%
%   lag              double (1 × (2*window+1))  
%       Time lag values in ms corresponding to each bin.
%
% Reference:
%   Adapted from the elegant phylib implementation here:
%   https://github.com/cortex-lab/phylib/blob/master/phylib/stats/ccg.py#L34
%
% Date:    20250704  
% Author:  Yue Huang

n_bins = floor(window/binwidth)+1;
auto_corr_right = zeros(1, n_bins); % the right side of auto_corr

shift = 1;
while true
    dt = spike_times(1+shift:end) - spike_times(1:end-shift);
    i_bin = int64(dt / binwidth) + 1;
    i_bin = i_bin(i_bin <= n_bins);

    if isempty(i_bin)
        break
    end

    counts = accumarray(i_bin(:), ones(length(i_bin),1))';

    auto_corr_right(1:length(counts)) = auto_corr_right(1:length(counts)) + counts;

    shift = shift+1;
end

auto_corr = [flip(auto_corr_right(2:end)), auto_corr_right];
lag = -(n_bins-1)*binwidth:binwidth:(n_bins-1)*binwidth;

assert(length(auto_corr) == length(lag));
assert(sum(lag==0) > 0);

end