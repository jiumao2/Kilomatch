function loss = lossFunLinearCorrection(params, sessions_pairs, dy, depth, linear_scale, n_session)
% params: 1 x n_session-2
% session_pairs: 2 x n
% dy: 1 x n
% depth: 1 x n
% depth_scale: 1 x 1, default: 1/1000
% n_session: 1 x 1

% extract the linear term and the constant term
p_linear = [0, params(1:n_session-1)];
p_constant = [0, params(n_session:end)];

p = linear_scale*p_linear(sessions_pairs).*depth + p_constant(sessions_pairs);
loss = sum(((p(2,:) - p(1,:)) - dy).^2, 'all');
end