function [auto_corr, lag] = computeAutoCorr(spike_times, window, binwidth)
% COMPUTEAUTOCORR compute the autocorrelogram using spike times without normalization
%
% spike_times: 1xn double array in ms
% window: 1x1 double in ms (default: 300)
% binwidth: 1x1 double in ms (default: 1)
%
% refer to the elegant Python impletantation from phylib:
%   https://github.com/cortex-lab/phylib/blob/master/phylib/stats/ccg.py#L34
%
% Converted to MATLAB by Yue Huang on 20250704
%

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