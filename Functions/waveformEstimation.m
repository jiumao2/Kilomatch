function waveform_out = waveformEstimation(waveform_mean, location, channel_locations, location_new)
% WAVEFORMESTIMATION  Interpolate waveform mean to a shifted probe location using 2D kernel methods.
%
% This function performs Kriging‐style interpolation of the mean waveform across recording channels.
% It computes a spatial kernel matrix between the original channel positions and those displaced
% by the difference between location_new and location. The resulting interpolation matrix is
% applied to waveform_mean to estimate the waveform at the new location.
%
% Inputs:
%   waveform_mean       double (n_channel × n_sample)
%       Mean waveform values across channels and time samples for one unit.
%   location            double (1 × 2)
%       [x, y] coordinates of the original waveform peak on the probe.
%   channel_locations   double (n_channel × 2)
%       [x, y] coordinates of each recording channel on the probe.
%   location_new        double (1 × 2)
%       [x, y] target coordinates for interpolation.
%
% Outputs:
%   waveform_out        double (n_channel × n_sample)
%       Interpolated waveform mean at the new probe location.
%
% Date:    20250821
% Author:  Yue Huang

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