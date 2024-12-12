function waveform_out = waveformEstimation(waveform_mean, location, chanMap, location_new, x, y, algorithm)
if nargin < 7
    algorithm = 'Krig';
end

% get n_nearest_channels from the channels with the largest peak-to-trough value
chanMap.xcoords = chanMap.xcoords(chanMap.connected == 1);
chanMap.ycoords = chanMap.ycoords(chanMap.connected == 1);

channel_locations = [chanMap.xcoords, chanMap.ycoords];

location_mapped_to_old = [x,y,0] - (location_new - location);

if strcmpi(algorithm, 'Raw')
    distance_to_location = sum((channel_locations - [x,y]).^2, 2);
    [~, idx_min] = min(distance_to_location);
    waveform_out = waveform_mean(idx_min,:);

    return
end

if strcmpi(algorithm, 'Nearest')
    distance_to_new_location = sum((channel_locations - location_mapped_to_old(1:2)).^2, 2);
    [~, idx_min] = min(distance_to_new_location);
    waveform_out = waveform_mean(idx_min,:);

    return
end

if strcmpi(algorithm, 'IDW')
    n_nearest_channels = 3;

    distance_to_new_location = sum((channel_locations - location_mapped_to_old(1:2)).^2, 2);
    [distance_this, idx_sort] = sort(distance_to_new_location);
    waveform_this = waveform_mean(idx_sort(1:n_nearest_channels),:);
    weight = 1./distance_this(1:n_nearest_channels);
    weight = weight./sum(weight);

    waveform_out = mean(waveform_this.*weight);
    return
end


if strcmpi(algorithm, 'Transformation')
    n_nearest_channels = 1;
    
    distance_to_location = sum((channel_locations - [x,y]).^2, 2);
    
    [~, idx_sorted] = sort(distance_to_location, 'ascend');
    idx_included = idx_sorted(1:n_nearest_channels);
    
    channel_locations_included = [channel_locations(idx_included, :), zeros(n_nearest_channels, 1)];
    
    weight0 = 1./sqrt(sum((channel_locations_included - location).^2, 2));
    weight1 = 1./sqrt(sum((location_mapped_to_old - location).^2));
    
    waveform_mean_included = waveform_mean(idx_included, :);
    
    waveform_out = mean(waveform_mean_included ./ weight0 .* weight1, 1);

    return
end

if strcmpi(algorithm, 'Krig')
    n_channels = 32;
    distance_to_location = sum((channel_locations - [x,y]).^2, 2);
    
    [~, idx_sorted] = sort(distance_to_location, 'ascend');
    idx_included = idx_sorted(1:n_channels);

    % 2D coordinates for interpolation 
    xp = channel_locations(idx_included,:);
    
    % 2D kernel of the original channel positions 
    Kxx = computeKernel2D(xp, xp);

    % 2D kernel of the new channel positions
%     yp = xp;
%     yp = yp - (location_new - location); % * sig;
%     Kyx = kernel2D(yp, xp);
    yp = location_mapped_to_old;
    Kyx = computeKernel2D(yp, xp);
    
    % kernel prediction matrix
    M = Kyx /(Kxx + .01 * eye(size(Kxx,1)));

    waveform_out = sum(waveform_mean(idx_included,:) .* M');

    return
end

end