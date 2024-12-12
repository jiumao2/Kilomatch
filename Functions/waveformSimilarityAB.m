function similarity = waveformSimilarityAB(spikeInfoA, spikeInfoB, n_channels, algorithm)
locationA = spikeInfoA.Location;
locationB = spikeInfoB.Location;

chanMap.xcoords = spikeInfoA.Xcoords;
chanMap.ycoords = spikeInfoA.Ycoords;
chanMap.kcoords = spikeInfoA.Kcoords;
chanMap.connected = ones(length(spikeInfoA.Xcoords), 1);

channel_locations = [chanMap.xcoords, chanMap.ycoords];
distance_to_locationB = sum((channel_locations - locationB(1:2)).^2, 2);

[~, idx_sorted] = sort(distance_to_locationB, 'ascend');
idx_included = idx_sorted(1:n_channels);

waveformsA = reshape(spikeInfoA.Waveform(idx_included, :)', 1, []);

waveforms_estimated = zeros(n_channels, size(spikeInfoA.Waveform, 2));
for k = 1:n_channels
    waveforms_estimated(k,:) = waveformEstimation(spikeInfoB.Waveform,...
        locationB,...
        chanMap,...
        locationA,...
        chanMap.xcoords(idx_included(k)), chanMap.ycoords(idx_included(k)),...
        algorithm);
end

waveformsB = reshape(waveforms_estimated', 1, []);

cc = corrcoef(waveformsA, waveformsB);
similarity = atanh(cc(1, 2));

end