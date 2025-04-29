function waveform_out = waveformEstimation(waveform_mean, location, channel_locations, location_new)
% Krigging interpolation on the reference probe

% get n_nearest_channels from the channels with the largest peak-to-trough value
location_mapped_to_old = channel_locations - (location_new(1:2) - location(1:2));

% 2D coordinates for interpolation 
xp = channel_locations;

% 2D kernel of the original channel positions 
Kxx = computeKernel2D(xp, xp);

% 2D kernel of the new channel positions
% yp = xp;
% yp = yp - (location_new - location); % * sig;
% Kyx = kernel2D(yp, xp);
yp = location_mapped_to_old;
Kyx = computeKernel2D(yp, xp);

% kernel prediction matrix
M = Kyx /(Kxx + .01 * eye(size(Kxx,1)));

waveform_out = M * waveform_mean;
end