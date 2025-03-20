function similarity = computeSimilarity(x, y)
    temp = corrcoef(x, y);
    similarity = atanh(temp(1,2));
    
    if isnan(similarity)
        similarity = 0;
    end
end