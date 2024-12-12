function similarity = PETH_Similarity(pethA, pethB)
temp = corrcoef(pethA, pethB);
similarity = atanh(temp(1,2));

if isnan(similarity)
    similarity = 0;
end

end