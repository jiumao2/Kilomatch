function similarity = waveformSimilarity(waveforms, waveform_channels, n_channels)
if nargin < 3
    n_channels = 36;
end

% similarityAB
channels = waveform_channels(1, 1:n_channels);
idx_channels_A = 1:n_channels;
idx_channels_B = zeros(1, n_channels);
for k = 1:n_channels
    temp = find(waveform_channels(2,:) == channels(k));
    if isempty(temp)
        similarity = 0;
        return
    end

    idx_channels_B(k) = temp;
end

waveformsA = reshape(waveforms(1, idx_channels_A, :), 1, []);
waveformsB = reshape(waveforms(2, idx_channels_B, :), 1, []);

temp = corrcoef(waveformsA, waveformsB);
similarityAB = temp(1,2);
similarityAB = atanh(similarityAB);

% similarityBA
channels = waveform_channels(2, 1:n_channels);
idx_channels_A = zeros(1, n_channels);
idx_channels_B = 1:n_channels;
for k = 1:n_channels
    temp = find(waveform_channels(1,:) == channels(k));
    if isempty(temp)
        similarity = 0;
        return
    end

    idx_channels_A(k) = temp;
end

waveformsA = reshape(waveforms(1, idx_channels_A, :), 1, []);
waveformsB = reshape(waveforms(2, idx_channels_B, :), 1, []);

temp = corrcoef(waveformsA, waveformsB);
similarityBA = temp(1,2);
similarityBA = atanh(similarityBA);

similarity = max(similarityAB, similarityBA);

end