function [nSame, nA, nB] = graphEditNumber(matA, matB)
% GRAPHEDITNUMBER  Compute number of "minimal edges" in two graphs and their overlap.
%
% Calculates edges in graph A and B and counts edges present in both.
% 1. Converts adjacency matrices matA and matB to MATLAB graph objects.
% 2. Finds connected components and sums component sizes minus one.
% 3. For the intersection graph (matA & matB), counts common edges.
%
% Inputs:
%   matA        logical matrix (n × n)
%       Adjacency matrix of graph A (symmetric, no self–loops).
%   matB        logical matrix (n × n), optional
%       Adjacency matrix of graph B; if omitted, matB = matA.
%
% Outputs:
%   nSame       integer scalar
%       Number of edges common to both graphs (A ∩ B).
%   nA          integer scalar
%       Total number of edges in graph A.
%   nB          integer scalar
%       Total number of edges in graph B.
%
% Date:    20250821  
% Author:  Yue Huang 

if nargin < 2
    matB = matA;
end

GA = graph(matA);
GB = graph(matB);
GAB = graph(matA & matB);

comp_A = conncomp(GA);
nA = sum(arrayfun(@(x)sum(comp_A == x)-1, 1:max(comp_A)));

comp_B = conncomp(GB);
nB = sum(arrayfun(@(x)sum(comp_B == x)-1, 1:max(comp_B)));

comp_AB = conncomp(GAB);
nSame = sum(arrayfun(@(x)sum(comp_AB == x)-1, 1:max(comp_AB)));

end