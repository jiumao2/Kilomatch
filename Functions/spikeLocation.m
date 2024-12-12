function [x, y, z, ptt] = spikeLocation(waveforms_mean, chanMap, n_nearest_channels, algorithm)
% waveforms_mean: n_channel x n_sample double
% chanMap: struct
% n_nearest_channels: 1 x 1 double about how many channels to include
% algorithm: 'center_of_mass' or 'monopolar_triangulation'
%
% monopolar_triangulation: refer to Boussard, Julien, Erdem Varol, Hyun Dong Lee, Nishchal Dethe, and Liam Paninski. “Three-Dimensional Spike Localization and Improved Motion Correction for Neuropixels Recordings.” In Advances in Neural Information Processing Systems, 34:22095–105. Curran Associates, Inc., 2021. https://proceedings.neurips.cc/paper/2021/hash/b950ea26ca12daae142bd74dba4427c8-Abstract.html.
% > https://spikeinterface.readthedocs.io/en/stable/modules/postprocessing.html#spike-locations

if nargin < 3
    n_nearest_channels = 20;
end
if nargin < 4
    algorithm = 'monopolar_triangulation';
end

% get n_nearest_channels from the channels with the largest peak-to-trough value
chanMap.xcoords = chanMap.xcoords(chanMap.connected == 1);
chanMap.ycoords = chanMap.ycoords(chanMap.connected == 1);

channel_locations = [chanMap.xcoords, chanMap.ycoords];

peaks_to_trough = max(waveforms_mean, [], 2) - min(waveforms_mean, [], 2);
[~, idx_max] = max(peaks_to_trough);

loc_max = channel_locations(idx_max, :);
distance_to_max = sum((channel_locations - loc_max).^2, 2);

[~, idx_sorted] = sort(distance_to_max, 'ascend');
idx_included = idx_sorted(1:n_nearest_channels);

% calculate the center_to_mass location
ptt_max = peaks_to_trough(idx_max);
ptt_this = peaks_to_trough(idx_included);
loc_this = channel_locations(idx_included,:);

loc_center_to_mass = sum(loc_this.*ptt_this, 1)./sum(ptt_this);

if strcmpi(algorithm, 'center_of_mass')
    x = loc_center_to_mass(1);
    y = loc_center_to_mass(2);
    z = 0;
    ptt = ptt_max;

    return
end

% calculate the monopolar_triangulation location

% % fminunc
% fun = @(x) sum(...
%     (ptt_this - x(4)./sqrt((loc_this(:,1)-x(1)).^2 + (loc_this(:,2)-x(2)).^2 + x(3).^2)).^2);
% 
% x0 = [loc_center_to_mass, 1, ptt_max];
% 
% options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e4);
% loc_monopolar_triangulation = fminunc(fun, x0, options);
% x = loc_monopolar_triangulation(1);
% y = loc_monopolar_triangulation(2);

% nonlinear least-square fitting
fun = @(x, loc_this) x(4)./sqrt((loc_this(:,1)-x(1)).^2 + (loc_this(:,2)-x(2)).^2 + x(3).^2);
x_bound_lower = [loc_center_to_mass(1)-1000, loc_center_to_mass(2)-1000, 0, 0];
x_bound_upper = [loc_center_to_mass(1)+1000, loc_center_to_mass(2)+1000, 1000*10, 1000*ptt_max];
x0 = [loc_center_to_mass, 1, ptt_max];

% disp('Calculating the location...');
options = optimoptions('lsqcurvefit', 'Display', 'off');
loc_monopolar_triangulation = lsqcurvefit(fun, x0, loc_this, ptt_this, x_bound_lower, x_bound_upper, options);
x = loc_monopolar_triangulation(1);
y = loc_monopolar_triangulation(2);
z = loc_monopolar_triangulation(3);
ptt = loc_monopolar_triangulation(4);
end