function PCs = getPCs(spikeInfo, nPCs)
if nargin < 2
    nPCs = 6;
end

all_waveforms = zeros(length(spikeInfo), size(spikeInfo(1).Waveform, 2));

for k = 1:length(spikeInfo)
    all_waveforms(k,:) = spikeInfo(k).Waveform(spikeInfo(k).Channel,:);
end

coeff = pca(all_waveforms);

PCs = coeff(:, 1:nPCs)';

end