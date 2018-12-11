function freq = sd_downsample_freq(cfg,freq)

% ----------------------------------------------------------------------- %
% This function downsamples Fieldtrip-structured time-frequency data into
% averaged time-bins.
% 
% DEFINE PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ----------------------------------------------------------------------- %

% check parameter
if isfield(freq,'powspctrm')
    datatype = 'freq';
elseif isfield(freq,'individual')
    datatype = 'timelock';
end

switch datatype
    
    case 'freq'

        % reduce time resolution to limit multiple comparisons
        nBins   = floor(round(diff(cfg.toi) ./ cfg.win_dur,3));
        tBins   = repmat([cfg.toi(1) cfg.toi(1)],[nBins,1]) + repmat([cfg.win_dur cfg.win_dur],[nBins,1]).*[0:nBins-1;1:nBins]';

        pow = nan(size(freq.powspctrm,1),size(freq.powspctrm,2),size(freq.powspctrm,3),nBins);
        for n = 1 : nBins
            tI = freq.time >= tBins(n,1) & freq.time <= tBins(n,2);
            pow(:,:,:,n) = nanmean(freq.powspctrm(:,:,:,tI),4);
        end

        freq.powspctrm = pow;
        freq.time = mean(tBins,2);
        
    case 'timelock'
        
        % reduce time resolution to limit multiple comparisons
        nBins   = floor(round(diff(cfg.toi) ./ cfg.win_dur,3));
        tBins   = repmat([cfg.toi(1) cfg.toi(1)],[nBins,1]) + repmat([cfg.win_dur cfg.win_dur],[nBins,1]).*[0:nBins-1;1:nBins]';

        pow = nan(size(freq.individual,1),size(freq.individual,2),nBins);
        for n = 1 : nBins
            tI = freq.time >= tBins(n,1) & freq.time <= tBins(n,2);
            pow(:,:,n) = nanmean(freq.individual(:,:,tI),3);
        end

        freq.individual = pow;
        freq.time = mean(tBins,2);
        
end