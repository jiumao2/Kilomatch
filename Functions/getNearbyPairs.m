function [idx_unit_pairs, session_pairs] = getNearbyPairs(max_distance, sessions, locations, Motion)
% GETNEARBYPairs  Find all pairs of units whose (corrected) depth difference  
%                 does not exceed a specified threshold.
%
% Identifies unit pairs with depth separation ≤ max_distance.  
% 1. Applies optional motion correction to raw depths.  
% 2. Computes pairwise absolute differences in corrected depth.  
% 3. Returns index pairs (i<j) and their corresponding session IDs.
%
% Inputs:
%   max_distance : scalar double
%       Maximum allowed depth difference in µm.
%   sessions     : integer vector (n_unit × 1)
%       Session index for each unit.
%   locations    : double matrix (n_unit × 2)
%       X and Y coordinates of each unit (Y = raw depth in µm).
%   Motion       : struct (optional)
%       Motion correction parameters with fields:
%         .LinearScale   scalar
%         .Linear        vector of length n_session
%         .Constant      vector of length n_session
%
% Outputs:
%   idx_unit_pairs : integer matrix (n_pairs × 2)
%       Each row [i, j], i<j, of unit indices within the depth threshold.
%   session_pairs  : integer matrix (n_pairs × 2)
%       Corresponding session indices [sessions(i), sessions(j)].
%
% Date:    20250821  
% Author:  Yue Huang  

if nargin < 4
    Motion = [];
end

n_unit = length(sessions);

if isempty(Motion) || isempty(Motion.Constant)
    corrected_locations = locations(:,2)';
else
    corrected_locations = zeros(1, n_unit);
    for k = 1:length(corrected_locations)
        corrected_locations(k) = locations(k,2)...
            - (Motion.LinearScale*Motion.Linear(sessions(k))*locations(k,2) + Motion.Constant(sessions(k)));
    end
end

y_distance_matrix = abs(corrected_locations - corrected_locations');

[idx_col, idx_row] = ind2sub(size(y_distance_matrix), 1:numel(y_distance_matrix));
idx_good = find(y_distance_matrix(:)' <= max_distance & idx_col > idx_row);
idx_unit_pairs = [idx_row(idx_good); idx_col(idx_good)]';

session_pairs = sessions(idx_unit_pairs);

end