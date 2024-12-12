function similarityAB = PC_SimilarityAB(spikeInfoA, spikeInfoB, PCs, n_channels, n_PC, algorithm)

location = spikeInfoA.Location(1:2);
channel_locations = [spikeInfoA.Xcoords, spikeInfoA.Ycoords];
distance = sum((channel_locations-location).^2, 2);

[~, idx_sort] = sort(distance);
idx_included = idx_sort(1:n_channels);

% chanMap
chanMap.xcoords = spikeInfoA.Xcoords;
chanMap.ycoords = spikeInfoA.Ycoords;
chanMap.kcoords = spikeInfoA.Kcoords;
chanMap.connected = ones(1, length(spikeInfoA.Xcoords));

PC_featuresA = zeros(n_channels, n_PC);
PC_featuresB = zeros(n_channels, n_PC);
for k = 1:length(idx_included)
    x = channel_locations(idx_included(k), 1);
    y = channel_locations(idx_included(k), 2);
    
    if strcmpi(algorithm, 'raw')
        waveform_corrected = spikeInfoB.Waveform(idx_included(k),:);
    else
        waveform_corrected = waveformEstimation(...
            spikeInfoB.Waveform,...
            spikeInfoB.Location,...
            chanMap,...
            spikeInfoA.Location,...
            x, y, algorithm);
    end

    for j = 1:n_PC
        PC_featuresA(k,j) = sum(spikeInfoA.Waveform(idx_included(k),:) .* PCs(j,:));
        PC_featuresB(k,j) = sum(waveform_corrected .* PCs(j,:));
    end
end

% temp = corrcoef(PC_featuresA(:), PC_featuresB(:));
% similarityAB = temp(1,2);
similarityAB = dot(PC_featuresA(:), PC_featuresB(:))./norm(PC_featuresA(:))./norm(PC_featuresB(:));
similarityAB = atanh(similarityAB);

end