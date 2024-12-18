function similarity = PETH_Similarity(spikeInfoA, spikeInfoB)
if ~isfield(spikeInfoA, 'PETH')
    similarity = NaN;
    return
end

temp = corrcoef(spikeInfoA.PETH, spikeInfoB.PETH);
similarity = atanh(temp(1,2));

if isnan(similarity)
    similarity = 0;
end

end