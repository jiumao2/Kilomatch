function loss = lossFunLinearCorrection(params, sessions_pairs, dy, depth, linear_scale, n_session)
% LOSSFUNLINEARCORRECTION  Loss function for linear session‐to‐session depth correction.
%
% Computes the sum of squared errors between predicted depth displacements and observed shifts.
% 1. Parses params into per‐session linear and constant offsets (session 1 fixed at zero).  
% 2. Predicts displacement p = linear_scale⋅linear(session)⋅depth_raw + constant(session).  
% 3. Calculates loss = Σ((p(2,:)–p(1,:) – dy).^2) across all unit pairs.
%
% Inputs:
%   params           : double vector (1 × [2*n_session – 2])  
%                      Concatenated linear (entries 2…n_session) and constant (entries n_session+1…2*n_session–2) terms; session 1 terms = 0.
%   sessions_pairs   : integer matrix (2 × n_pairs)  
%                      Session indices [sess_i; sess_j] for each observed pair.
%   dy               : double vector (1 × n_pairs)  
%                      Observed depth difference (depth_j – depth_i) for each pair.
%   depth            : double vector (1 × n_pairs)  
%                      Raw depth in µm of the first‐session unit in each pair.
%   linear_scale     : double scalar  
%                      Global scaling factor applied to linear terms.
%   n_session        : integer scalar  
%                      Total number of recording sessions.
%
% Outputs:
%   loss             : double scalar  
%                      Sum of squared errors between predicted and observed displacements.
%
% Date:    20250821  
% Author:  Yue Huang 

% extract the linear term and the constant term
p_linear = [0, params(1:n_session-1)];
p_constant = [0, params(n_session:end)];

p = linear_scale*p_linear(sessions_pairs).*depth + p_constant(sessions_pairs);
loss = sum(((p(2,:) - p(1,:)) - dy).^2, 'all');
end