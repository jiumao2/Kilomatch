function loss = lossFunLinearDefault(params, sessions_pairs, dy)
% params: 1 x n_session
% session_pairs: 2 x n
% dy: 1 x n

loss = sum(((params(sessions_pairs(2,:)) - params(sessions_pairs(1,:))) - dy).^2, 'all');
end