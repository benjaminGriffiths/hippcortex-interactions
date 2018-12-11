function freq = sd_ztransform_freq(freq)

% -------------------------------------------------------------- %
% this function z-transforms power across the entire time series %
% -------------------------------------------------------------- %

avgPow = repmat(nanmean(nanmean(freq.powspctrm,4),1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
stdPow = repmat(nanstd(nanmean(freq.powspctrm,4),[],1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
freq.powspctrm = (freq.powspctrm - avgPow) ./ stdPow;