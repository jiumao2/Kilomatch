function similarity = ISI_Similarity(spikeInfoA, spikeInfoB)
cc = corrcoef(spikeInfoA.ISI, spikeInfoB.ISI);
similarity = atanh(cc(1, 2)); 

if isnan(similarity)
    similarity = 0;
end

end