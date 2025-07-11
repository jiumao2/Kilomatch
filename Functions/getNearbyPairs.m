function [idx_unit_pairs, session_pairs] = getNearbyPairs(max_distance, sessions, locations, Motion)
% max_distance: 1x1 scalar, um
% sessions: n_unit x 1 array
% locations: n_unit x 2 array
% motion: 1 x n_session array

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