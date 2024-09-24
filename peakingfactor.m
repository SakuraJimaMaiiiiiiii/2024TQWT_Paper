function peakingfactor=peakingfactor(x)
    peakingfactor=0.5*(max(x)-min(x))/rms(x);
end

