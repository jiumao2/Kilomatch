function loss = lossFunLinearDefault(params, sessions_pairs, dy)
% LOSSFUNLINEARDEFAULT  Loss function for default linear session depth correction.
%
% Computes squared‐error loss between session‐to‐session parameter offsets
% and observed depth differences.
%
% Inputs:
%   params             double vector (1 × n_session)
%       Per‐session depth offset parameters (session 1…n_session).
%   sessions_pairs     integer matrix (2 × n_pairs)
%       Session index pairs [sess_i; sess_j] for each observation.
%   dy                 double vector (1 × n_pairs)
%       Observed depth difference (depth_j – depth_i) for each pair.
%
% Outputs:
%   loss               double scalar
%       Sum of squared errors between predicted and observed depth shifts.
%
% Date:    20250821  
% Author:  Yue Huang  

loss = sum(((params(sessions_pairs(2,:)) - params(sessions_pairs(1,:))) - dy).^2, 'all');

end