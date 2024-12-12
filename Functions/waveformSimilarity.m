function similarity = waveformSimilarity(spikeInfoA, spikeInfoB, n_channels, algorithm)
if nargin < 3
    n_channels = 32;
end

if nargin < 4
    algorithm = 'Krig';
end

similarity = max(...
    waveformSimilarityAB(spikeInfoA, spikeInfoB, n_channels, algorithm),...
    waveformSimilarityAB(spikeInfoB, spikeInfoA, n_channels, algorithm));

end



