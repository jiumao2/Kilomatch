function [nSame, nA, nB] = graphEditNumber(matA, matB)
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