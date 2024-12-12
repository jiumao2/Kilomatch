function similarity = PC_Similarity(spikeInfoA, spikeInfoB, PCs, n_channels, n_PC, algorithm)
if nargin < 4
    n_channels = 36;
end

if nargin < 5
    n_PC = size(PCs, 1);
end

if nargin < 6
    algorithm = 'Krig';
end

similarity = max(...
    PC_SimilarityAB(spikeInfoA, spikeInfoB, PCs, n_channels, n_PC, algorithm),...
    PC_SimilarityAB(spikeInfoB, spikeInfoA, PCs, n_channels, n_PC, algorithm));

end
