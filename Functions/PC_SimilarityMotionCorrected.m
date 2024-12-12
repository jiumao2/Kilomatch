function similarity = PC_SimilarityMotionCorrected(PC_features, PC_channels, n_channels, nPC)
if nargin < 3
    n_channels = 36;
end

if nargin < 4
    nPC = size(PC_features, 3);
end

% similarityAB
channels = PC_channels(1, 1:n_channels);
idx_channels_A = 1:n_channels;
idx_channels_B = zeros(1, n_channels);
for k = 1:n_channels
    temp = find(PC_channels(2,:) == channels(k));
    if isempty(temp)
        similarity = 0;
        return
    end

    idx_channels_B(k) = temp;
end

PC_featuresA = reshape(PC_features(1, idx_channels_A, 1:nPC), 1, []);
PC_featuresB = reshape(PC_features(2, idx_channels_B, 1:nPC), 1, []);
similarityAB = dot(PC_featuresA, PC_featuresB)./norm(PC_featuresA)./norm(PC_featuresB);
% temp = corrcoef(PC_featuresA, PC_featuresB);
% similarityAB = temp(1,2);
similarityAB = atanh(similarityAB);

% similarityBA
channels = PC_channels(2, 1:n_channels);
idx_channels_A = zeros(1, n_channels);
idx_channels_B = 1:n_channels;
for k = 1:n_channels
    temp = find(PC_channels(1,:) == channels(k));
    if isempty(temp)
        similarity = 0;
        return
    end

    idx_channels_A(k) = temp;
end

PC_featuresA = reshape(PC_features(1, idx_channels_A, 1:nPC), 1, []);
PC_featuresB = reshape(PC_features(2, idx_channels_B, 1:nPC), 1, []);
similarityBA = dot(PC_featuresA, PC_featuresB)./norm(PC_featuresA)./norm(PC_featuresB);
% temp = corrcoef(PC_featuresA, PC_featuresB);
% similarityBA = temp(1,2);
similarityBA = atanh(similarityBA);

similarity = max(similarityAB, similarityBA);

end