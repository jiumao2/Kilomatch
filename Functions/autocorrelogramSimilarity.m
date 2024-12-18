function similarity = autocorrelogramSimilarity(spikeInfoA, spikeInfoB)
cc = corrcoef(spikeInfoA.AutoCorr, spikeInfoB.AutoCorr);
similarity = atanh(cc(1,2));

if isnan(similarity)
    similarity = 0;
end

end