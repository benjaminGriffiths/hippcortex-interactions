%% Cross-Correlation Analysis
% This script aims to identify the relationship between neocortical (NC) 
% and hippocampal power during successful memory formation and retrieval.
clearvars
close all
clc

% define key directory
home_dir = 'E:\bjg335\projects\hippcortex-interactions';
data_dir = 'Y:\projects\intracranial_sync_desync';

% load contact labels
load([home_dir,'\contact_locations.mat'])

% load peak frequencies
load([home_dir,'\peak_frequencies.mat'])

% add subfunctions
addpath([home_dir,'\additional_functions'])

% define key parameters
subj_no     = [1 2 3 4 6 7 8];

%% Get Time-Frequency Spectrum
% cycle through every participant
for subj = subj_no
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);
    
    % pre-define cell for peak time-frequency data
    freq_peak = cell(3,numel(data));
    
    % predefine time window of interest based on encoding/retrieval data
    toi{1} = 3.5 : 0.01 : 4.5;
    toi{2} = 0.5 : 0.01 : 1.5;
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(data)
    
        % calculate time-frequency
        cfg             = [];
        cfg.keeptrials  = 'yes';
        cfg.method      = 'wavelet'; 
        cfg.width       = 5; 
        cfg.output      = 'pow';	
        cfg.pad         = 'nextpow2';
        cfg.foi         = 1:0.5:100;          
        cfg.toi         = toi{di};
        freq            = ft_freqanalysis(cfg, data{di});

        % subtract 1/f noise
        cfg             = [];
        cfg.toi         = [freq.time(1) freq.time(end)];
        freq            = sd_subtr1of(cfg, freq);

        % smooth data
        cfg             = [];
        cfg.fwhm_t      = 0.05;
        cfg.fwhm_f      = 1;
        freq            = smooth_TF_GA(cfg,freq);
        
        % convert power spectrum to single to reduce RAM load
        freq.powspctrm = single(freq.powspctrm);
 
        % define ATL alpha
        cfg                     = [];
        cfg.frequency           = peak_frequencies.alpha(subj,:);
        cfg.channel             = freq.label(contact_locations.atl{subj,:});
        cfg.avgoverfreq         = 'yes';
        freq_peak{1,di}         = ft_selectdata(cfg,freq);
        freq_peak{1,di}.freq    = 10;

        % define PTPR alpha
        cfg                     = [];
        cfg.frequency           = peak_frequencies.alpha(subj,:);
        cfg.channel             = freq.label(contact_locations.nc{subj,:});
        cfg.avgoverfreq         = 'yes';
        freq_peak{2,di}         = ft_selectdata(cfg,freq);
        freq_peak{2,di}.freq    = 10;

        % define hippocampal theta
        cfg                     = [];
        cfg.frequency           = peak_frequencies.theta(subj,:);
        cfg.channel             = freq.label(contact_locations.hippo{subj,:});
        cfg.avgoverfreq         = 'yes';
        tmp{1,1}                = ft_selectdata(cfg,freq);
        tmp{1,1}.freq           = 4;

        % define hippocampal encoding gamma
        cfg                     = [];
        cfg.frequency           = peak_frequencies.enc_gamma(subj,:);
        cfg.channel             = freq.label(contact_locations.hippo{subj,:});
        cfg.avgoverfreq         = 'yes';
        tmp{2,1}                = ft_selectdata(cfg,freq);
        tmp{2,1}.freq           = 60;

        % define hippocampal retrieval gamma
        cfg                     = [];
        cfg.frequency           = peak_frequencies.ret_gamma(subj,:);
        cfg.channel             = freq.label(contact_locations.hippo{subj,:});
        cfg.avgoverfreq         = 'yes';
        tmp{3,1}                = ft_selectdata(cfg,freq);
        tmp{3,1}.freq           = 45;
    
        % concatenate hippocampal data across frequencies
        cfg                     = [];
        cfg.parameter           = 'powspctrm';
        cfg.appenddim           = 'freq';
        freq_peak{3,di}         = ft_appendfreq(cfg, tmp{:});
    end
    
    % save data
    save([data_dir,'\data\xc\pp',num2str(subj),'_freq_peak.mat'],'freq_peak');
    
    % clear all non-essential variables
    keep contact_locations data_dir home_dir peak_frequencies subj_no
end
    
%% Calculate Cross-Correlation
% cycle through every participant
for subj = subj_no
    
    % load data
    load([data_dir,'\data\xc\pp',num2str(subj),'_freq_peak.mat'])
        
    % predefine frequency cell
    freq_xc = cell(size(freq_peak,2),1);
    
    % define cross-corelation parameters
    xc_lag = 30;   % lag in samples
    xc_tD  = 0.01; % duration of samples
    
    % cycle through encoding and retrieval data
    for di = 1 : size(freq_peak,2)
        
        % create cross-correlation data structure
        freq_xc{di} = struct('cfg',[],...
                         'freq',[4 60 45],...
                         'time',linspace(-xc_lag.*xc_tD,xc_lag.*xc_tD,xc_lag.*2+1),...
                         'label',{{'atl_hipp','ptl_hipp'}},...
                         'powspctrm',nan(size(freq_peak{1,di}.powspctrm,1),2,3,xc_lag.*2+1),...
                         'dimord','rpt_chan_freq_time',...
                         'trialinfo',freq_peak{1,di}.trialinfo);
        
        % cycle through NC ROI
        for ri = 1 : numel(freq_xc{di}.label)
            
            % determine ATL-Hippo. contact pairs
            chancmb = sd_get_pairs(freq_peak{ri,di},freq_peak{3,di});

            % pre-define matrix for cross-correlation
            XC_tmp = nan(size(freq_peak{ri,di}.powspctrm,1),size(chancmb,1),numel(freq_peak{ri,di}.freq),xc_lag.*2+1);

            % cycle through each trial
            for trial = 1 : size(freq_peak{ri,di}.powspctrm,1)

                % cycle through each channel combination
                for chan = 1 : size(chancmb,1)

                    % extract NC and hippocampal time-series
                    signal_A = squeeze(freq_peak{ri,di}.powspctrm(trial,chancmb(chan,1),:,:));
                    signal_B = squeeze(freq_peak{3,di}.powspctrm(trial,chancmb(chan,2),:,:));

                    % rotate to trls_freqs if necessary
                    if size(signal_B,1) < size(signal_B,2); signal_B = signal_B'; end
                    
                    % cycle through each hippocampal frequency
                    for freq = 1 : size(signal_B,2)

                        % calculate XC (using atanh to normalise Pearson's coefficent)
                        XC_tmp(trial,chan,freq,:) = atanh(crosscorr(signal_B(:,freq),signal_A,xc_lag));

                    end
                end
            end

            % average XC values over contact combinations and add to structure
            freq_xc{di}.powspctrm(:,ri,:,:) = mean(XC_tmp,2);

        end
    end
    
    % update command line
    fprintf('Subject %d cross-correlation complete...\n',subj)
    
    % save data
    save([data_dir,'\data\xc\pp',num2str(subj),'_freq_xc.mat'],'freq_xc')
    
    % clear all non-essential variables
    keep contact_locations data_dir home_dir peak_frequencies subj_no
end

%% Create Memory Contrast
% cycle through every participant
for subj = subj_no
    
    % load data
    load([data_dir,'\data\xc\pp',num2str(subj),'_freq_xc.mat'])
     
    % predefine cell for freq data
    freq_hits = cell(2,1);
    freq_misses = cell(2,1);
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(freq_xc)
        
        % split trials based on memory performance
        cfg             = [];
        cfg.trials      = freq_xc{di}.trialinfo == 1;
        freq_hits{di}   = ft_selectdata(cfg,freq_xc{di});

        cfg             = [];
        cfg.trials      = freq_xc{di}.trialinfo == 0;
        freq_misses{di} = ft_selectdata(cfg,freq_xc{di});        

        % average over trials
        cfg             = [];
        freq_hits{di}   = ft_freqdescriptives(cfg, freq_hits{di});
        freq_misses{di} = ft_freqdescriptives(cfg, freq_misses{di});
    end
     
    % save data
    save([data_dir,'\data\xc\pp',num2str(subj),'_freq_sme.mat'],'freq_hits','freq_misses')
    
    % clear all non-essential variables
    keep contact_locations data_dir home_dir peak_frequencies subj_no
end

%% Calculate Group Averages
% create group cell
group_hits      = cell(max(subj_no),2);
group_misses    = cell(max(subj_no),2);

% cycle through each participant
for subj = subj_no
    
    % load data
    load([data_dir,'\data\xc\pp',num2str(subj),'_freq_sme.mat']);
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(freq_hits)
        group_hits{subj,di} = freq_hits{di};
        group_misses{subj,di} = freq_misses{di};
    end
end

% get grand average
grand_hits          = cell(2,1);
grand_misses        = cell(2,1);
cfg                 = [];
cfg.keepindividual  = 'yes';
for di = 1 : size(group_hits,2)
    grand_hits{di,1}    = ft_freqgrandaverage(cfg,group_hits{subj_no,di});
    grand_misses{di,1}  = ft_freqgrandaverage(cfg,group_misses{subj_no,di});
end

% save
save([data_dir,'\data\xc\grand_freq.mat'],'grand_hits','grand_misses')
      
% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Run Inferential Statistics
% load data
load([data_dir,'\data\xc\grand_freq.mat'])

% define frequency of interest
foi = [59 61; 44 46];

% predefine matrices to hold statisitcal data
p = zeros(2,6);
d = zeros(2,6);

% cycle through both data types
for di = 1 : numel(grand_hits)
    
    % downsample time domain
    cfg                 = [];
    cfg.win_dur         = 0.1;
    cfg.toi             = [grand_hits{di}.time(1) grand_hits{di}.time(end)];
    grand_hits{di}      = sd_downsample_freq(cfg,grand_hits{di});
    grand_misses{di}    = sd_downsample_freq(cfg,grand_misses{di});
    
    % create design matrix
    design      = [];
    design(1,:) = [1:size(grand_hits{di}.powspctrm,1), 1:size(grand_hits{di}.powspctrm,1)];
    design(2,:) = [ones(1,size(grand_hits{di}.powspctrm,1)) ones(1,size(grand_hits{di}.powspctrm,1))+1];

    % define statistic configuration
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'ft_statfun_depsamplesT';
    cfg.correctm            = 'no';
    cfg.alpha               = 0.05;
    cfg.numrandomization    = 5000;
    cfg.design              = design;
    cfg.ivar                = 2;
    cfg.uvar                = 1;
    cfg.tail                = -1;
       
    % --- test alpha/gamma correlation --- %
    cfg.channel             = 'atl_hipp';
    cfg.frequency           = foi(di,:);
    stat                    = ft_freqstatistics(cfg, grand_hits{di}, grand_misses{di});

    % extract p-value
    [~,~,p(di,:)] = fdr(squeeze(stat.prob));
    
    % extract cohen's d
    for t = 1 : size(stat.stat,3)
        d(di,t) = computeCohen_d(squeeze(grand_hits{di}.powspctrm(:,1,3,t)),...
                              squeeze(grand_misses{di}.powspctrm(:,1,3,t)),'paired');
    end   
    
    % calculate number of electrodes showing effect
    if di == 1
        percent_trending(di,1) = sum((grand_hits{di}.powspctrm(:,1,3,2) - grand_misses{di}.powspctrm(:,1,3,2))<0) ./ size(grand_hits{di}.powspctrm,1);
    else
        percent_trending(di,1) = sum((grand_hits{di}.powspctrm(:,1,2,6) - grand_misses{di}.powspctrm(:,1,2,6))<0) ./ size(grand_hits{di}.powspctrm,1);
    end
end

% save
save([data_dir,'\data\xc\stat_hitMissContrast.mat'],'p','d','percent_trending')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Contrast 60Hz Encoding with 45Hz Retrieval
% load data
load([data_dir,'\data\xc\grand_freq.mat'])

% cycle through both data types
for di = 1 : numel(grand_hits)
    
    % get difference between hits and misses
    grand_diff{di} = grand_hits{di}; %#ok<SAGROW>
    grand_diff{di}.powspctrm = grand_hits{di}.powspctrm - grand_misses{di}.powspctrm;%#ok<SAGROW>
    
    % downsample time domain
    cfg                 = [];
    cfg.win_dur         = 0.1;
    cfg.toi             = [grand_diff{di}.time(1) grand_diff{di}.time(end)];
    grand_diff{di}      = sd_downsample_freq(cfg,grand_diff{di}); %#ok<SAGROW>      
end

% select relevant frequency
cfg             = [];
cfg.frequency   = [60 60];
grand_freq{1}   = ft_selectdata(cfg,grand_diff{1});

cfg             = [];
cfg.frequency   = [45 45];
grand_freq{2}   = ft_selectdata(cfg,grand_diff{2});

% rename frequency for consistency
grand_freq{1}.freq = 50;
grand_freq{2}.freq = 50;

% create design matrix
design      = [];
design(1,:) = [1:size(grand_freq{di}.powspctrm,1), 1:size(grand_freq{di}.powspctrm,1)];
design(2,:) = [ones(1,size(grand_freq{di}.powspctrm,1)) ones(1,size(grand_freq{di}.powspctrm,1))+1];

% define statistic configuration
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'no';
cfg.alpha               = 0.05;
cfg.numrandomization    = 5000;
cfg.design              = design;
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.channel             = 'atl_hipp';

% run one-sample statistics
stat = ft_freqstatistics(cfg, grand_freq{1}, grand_freq{2});

% extract p-values
[~,~,p] = fdr(squeeze(stat.prob));

% extract cohen's d
for t = 1 : size(stat.stat,3)
    d(t) = computeCohen_d(squeeze(grand_freq{1}.powspctrm(:,1,1,t)),...
                          squeeze(grand_freq{2}.powspctrm(:,1,1,t)),'paired');
end    

% calculate number of electrodes showing effect
percent_trending(1,1) = sum((grand_freq{1}.powspctrm(:,1,1,2) - grand_freq{2}.powspctrm(:,1,1,2))<0) ./ size(grand_freq{di}.powspctrm,1);
percent_trending(2,1) = sum((grand_freq{1}.powspctrm(:,1,1,6) - grand_freq{2}.powspctrm(:,1,1,6))>0) ./ size(grand_freq{di}.powspctrm,1);

% save
save([data_dir,'\data\xc\stat_encRetContrast.mat'],'p','d','percent_trending')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Run Random Effects ANOVA
% predefine vectors for ANOVA analysis
X = []; % power
Y = []; % subject
A = []; % frequency
B = []; % encoding/retrieval

% load data
load([data_dir,'\data\xc\grand_freq.mat'])

% cycle through both data types
for di = 1 : numel(grand_hits)
    
    % downsample time domain
    cfg                 = [];
    cfg.win_dur         = 0.1;
    cfg.toi             = [grand_hits{di}.time(1) grand_hits{di}.time(end)];
    grand_hits{di}      = sd_downsample_freq(cfg,grand_hits{di});
    grand_misses{di}    = sd_downsample_freq(cfg,grand_misses{di});
    
    % get difference between hits and misses
    grand_diff{di} = grand_hits{di};
    grand_diff{di}.powspctrm = grand_hits{di}.powspctrm - grand_misses{di}.powspctrm;

end

% define encoding and retrieval time index
time_idx = [2 6];

% cycle through each participant
for subj = 1 : numel(subj_no)

    % cycle through both data types
    for di = 1 : numel(grand_hits)

        % extract power (ATL/Gamma/Pre-Post) and add to group vector
        X(end+1:end+2,1) = grand_diff{di}.powspctrm(subj,1,[2 3],time_idx(di));
    
        % add frequency for subject
        A(end+1:end+2,1) = [1 2];

        % add trial type for subject
        B(end+1:end+2,1) = [di di];

        % add subject number
        Y(end+1:end+2,1) = [subj subj];
    end    
end

% run 2x2 RM-ANOVA
p = rm_anova2_np(X,Y,A,B,5000);

% combined data into single matrix
E = [X A B];
[~,idx] = sort(A);
E = E(idx,:);

Ei = E(E(:,2) == 1,:);
[~,idx] = sort(Ei(:,3));
Ei = Ei(idx,:);

Ej = E(E(:,2) == 2,:);
[~,idx] = sort(Ej(:,3));
Ej = Ej(idx,:);

E = [Ei;Ej];

% get effect size
e = mes2way(E(:,1),E(:,2:3),'partialeta2');

% save
save([data_dir,'\data\xc\stat_gammaRE.mat'],'p','e')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations subj_no file_label

%% Prepare Data for Fixed Effects ANOVA
% predefine matrix for group data
group_xc{1,1} = [];
group_xc{1,2} = [];
group_xc{2,1} = [];
group_xc{2,2} = [];

% load contact labels
load([home_dir,'\contact_locations.mat'])

% cycle through every participant
for subj = subj_no
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);
    
    % pre-define cell for peak time-frequency data
    freq_peak = cell(3,numel(data));
    
    % predefine time window of interest based on encoding/retrieval data
    toi{1} = 3.5 : 0.01 : 4.5;
    toi{2} = 0.5 : 0.01 : 1.5;
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(data)
    
        % calculate time-frequency
        cfg             = [];
        cfg.keeptrials  = 'yes';
        cfg.method      = 'wavelet'; 
        cfg.width       = 5; 
        cfg.output      = 'pow';	
        cfg.pad         = 'nextpow2';
        cfg.foi         = 1:0.5:100;          
        cfg.toi         = toi{di};
        freq            = ft_freqanalysis(cfg, data{di});

        % subtract 1/f noise
        cfg             = [];
        cfg.toi         = [freq.time(1) freq.time(end)];
        freq            = sd_subtr1of(cfg, freq);

        % smooth data
        cfg             = [];
        cfg.fwhm_t      = 0.05;
        cfg.fwhm_f      = 1;
        freq            = smooth_TF_GA(cfg,freq);
        
        % convert power spectrum to single to reduce RAM load
        freq.powspctrm = single(freq.powspctrm);
 
        % define ATL alpha
        cfg                     = [];
        cfg.frequency           = peak_frequencies.alpha(subj,:);
        cfg.channel             = freq.label(contact_locations.atl{subj,:});
        cfg.avgoverfreq         = 'yes';
        freq_peak{1,di}         = ft_selectdata(cfg,freq);
        freq_peak{1,di}.freq    = 10;

        % define PTPR alpha
        cfg                     = [];
        cfg.frequency           = peak_frequencies.alpha(subj,:);
        cfg.channel             = freq.label(contact_locations.nc{subj,:});
        cfg.avgoverfreq         = 'yes';
        freq_peak{2,di}         = ft_selectdata(cfg,freq);
        freq_peak{2,di}.freq    = 10;

        % define hippocampal theta
        cfg                     = [];
        cfg.frequency           = peak_frequencies.theta(subj,:);
        cfg.channel             = freq.label(contact_locations.hippo{subj,:});
        cfg.avgoverfreq         = 'yes';
        tmp{1,1}                = ft_selectdata(cfg,freq);
        tmp{1,1}.freq           = 4;

        % define hippocampal encoding gamma
        cfg                     = [];
        cfg.frequency           = peak_frequencies.ret_gamma(subj,:);
        cfg.channel             = freq.label(contact_locations.hippo{subj,:});
        cfg.avgoverfreq         = 'yes';
        tmp{2,1}                = ft_selectdata(cfg,freq);
        tmp{2,1}.freq           = 45;

        % define hippocampal retrieval gamma
        cfg                     = [];
        cfg.frequency           = peak_frequencies.enc_gamma(subj,:);
        cfg.channel             = freq.label(contact_locations.hippo{subj,:});
        cfg.avgoverfreq         = 'yes';
        tmp{3,1}                = ft_selectdata(cfg,freq);
        tmp{3,1}.freq           = 60;
    
        % concatenate hippocampal data across frequencies
        cfg                     = [];
        cfg.parameter           = 'powspctrm';
        cfg.appenddim           = 'freq';
        freq_peak{3,di}         = ft_appendfreq(cfg, tmp{:});
    end
        
    % predefine frequency cell
    freq_xc = cell(size(freq_peak,2),1);
    
    % define cross-corelation parameters
    xc_lag = 30;   % lag in samples
    xc_tD  = 0.01; % duration of samples
    
    % cycle through encoding and retrieval data
    for di = 1 : size(freq_peak,2)
        
        % determine ATL-Hippo. contact pairs
        chancmb = sd_get_pairs(freq_peak{1,di},freq_peak{3,di});

        % pre-define matrix for cross-correlation
        XC_tmp = nan(size(freq_peak{1,di}.powspctrm,1),size(chancmb,1),numel(freq_peak{1,di}.freq),xc_lag.*2+1);

        % cycle through each trial
        for trial = 1 : size(freq_peak{1,di}.powspctrm,1)

            % cycle through each channel combination
            for chan = 1 : size(chancmb,1)

                % extract NC and hippocampal time-series
                signal_A = squeeze(freq_peak{1,di}.powspctrm(trial,chancmb(chan,1),:,:));
                signal_B = squeeze(freq_peak{3,di}.powspctrm(trial,chancmb(chan,2),:,:));

                % rotate to trls_freqs if necessary
                if size(signal_B,1) < size(signal_B,2); signal_B = signal_B'; end

                % cycle through each hippocampal frequency
                for freq = 1 : size(signal_B,2)

                    % calculate XC (using atanh to normalise Pearson's coefficent)
                    XC_tmp(trial,chan,freq,:) = atanh(crosscorr(signal_B(:,freq),signal_A,xc_lag));

                end
            end
        end
        
        % calculate hits and misses
        xc_hits 	= squeeze(mean(XC_tmp(freq_peak{1,di}.trialinfo==1,:,:,:)));
        xc_misses   = squeeze(mean(XC_tmp(freq_peak{1,di}.trialinfo==0,:,:,:)));
        
        % add to group data
        group_xc{di,1}(end+1:end+size(xc_hits,1),:,:) = xc_hits;
        group_xc{di,2}(end+1:end+size(xc_misses,1),:,:) = xc_misses;
        clear xc_diff
    end
    
    % update command line
    fprintf('Subject %d cross-correlation complete...\n',subj)
    
    % clear all non-essential variables
    keep contact_locations data_dir home_dir peak_frequencies subj_no group_xc
end

% create data structure
for i = 1 : 2
    for j = 1 : 2
        grand_FE{i,j} = struct('cfg',[],...
                             'freq',[4 45 60],...
                             'time',-0.3:0.01:0.3,...
                             'label',{{'dummy'}},...
                             'dimord','subj_chan_freq_time',...
                             'powspctrm',permute(group_xc{i,j},[1 4 2 3])); %#ok<SAGROW>
    end
end

% save
save([data_dir,'\data\xc\fixedEffectsPow.mat'],'grand_FE')

%% Run Time-Series Inferential Statistics
% load data
load([data_dir,'\data\xc\fixedEffectsPow.mat'])

% define frequency of interest
foi = [59 61; 44 46];

% predefine matrices to hold statisitcal data
p = zeros(2,6);
d = zeros(2,6);
percent_trending = zeros(2,1);

% cycle through both data types
for di = 1 : size(grand_FE,1)
    
    % downsample time domain
    cfg                 = [];
    cfg.win_dur         = 0.1;
    cfg.toi             = [grand_FE{di,1}.time(1) grand_FE{di,1}.time(end)];
    grand_hits{di}      = sd_downsample_freq(cfg,grand_FE{di,1});
    grand_misses{di}    = sd_downsample_freq(cfg,grand_FE{di,2});
    
    % create design matrix
    design      = [];
    design(1,:) = [1:size(grand_hits{di}.powspctrm,1), 1:size(grand_hits{di}.powspctrm,1)];
    design(2,:) = [ones(1,size(grand_hits{di}.powspctrm,1)) ones(1,size(grand_hits{di}.powspctrm,1))+1];

    % define statistic configuration
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'ft_statfun_depsamplesT';
    cfg.correctm            = 'no';
    cfg.alpha               = 0.05;
    cfg.numrandomization    = 5000;
    cfg.design              = design;
    cfg.ivar                = 2;
    cfg.uvar                = 1;
    cfg.tail                = -1;
       
    % --- test alpha/gamma correlation --- %
    cfg.frequency           = foi(di,:);
    stat                    = ft_freqstatistics(cfg, grand_hits{di}, grand_misses{di});

    % extract p-value
    [~,~,p(di,:)] = fdr(squeeze(stat.prob));
    
    % extract cohen's d
    for t = 1 : size(stat.stat,3)
        d(di,t) = computeCohen_d(squeeze(grand_hits{di}.powspctrm(:,1,3,t)),...
                              squeeze(grand_misses{di}.powspctrm(:,1,3,t)),'paired');
    end    
    
    % calculate number of electrodes showing effect
    if di == 1
        percent_trending(di,1) = sum((grand_hits{di}.powspctrm(:,:,3,2) - grand_misses{di}.powspctrm(:,:,3,2))<0) ./ size(grand_hits{di}.powspctrm,1);
    else
        percent_trending(di,1) = sum((grand_hits{di}.powspctrm(:,:,2,6) - grand_misses{di}.powspctrm(:,:,2,6))<0) ./ size(grand_hits{di}.powspctrm,1);
    end
end

% save
save([data_dir,'\data\xc\stat_hitMissFEContrast.mat'],'p','d','percent_trending')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Contrast 60Hz Encoding with 45Hz Retrieval
% load data
load([data_dir,'\data\xc\fixedEffectsPow.mat'])

% cycle through both data types
for di = 1 : size(grand_FE,1)
    
    % get difference between hits and misses
    grand_diff{di} = grand_FE{di,1}; %#ok<SAGROW>
    grand_diff{di}.powspctrm = grand_FE{di,1}.powspctrm - grand_FE{di,2}.powspctrm;%#ok<SAGROW>
    
    % downsample time domain
    cfg                 = [];
    cfg.win_dur         = 0.1;
    cfg.toi             = [grand_diff{di}.time(1) grand_diff{di}.time(end)];
    grand_diff{di}      = sd_downsample_freq(cfg,grand_diff{di}); %#ok<SAGROW>      
end

% select relevant frequency
cfg             = [];
cfg.frequency   = [60 60];
grand_freq{1}   = ft_selectdata(cfg,grand_diff{1});

cfg             = [];
cfg.frequency   = [45 45];
grand_freq{2}   = ft_selectdata(cfg,grand_diff{2});

% rename frequency for consistency
grand_freq{1}.freq = 50;
grand_freq{2}.freq = 50;

% create design matrix
design      = [];
design(1,:) = [1:size(grand_freq{di}.powspctrm,1), 1:size(grand_freq{di}.powspctrm,1)];
design(2,:) = [ones(1,size(grand_freq{di}.powspctrm,1)) ones(1,size(grand_freq{di}.powspctrm,1))+1];

% define statistic configuration
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'no';
cfg.alpha               = 0.05;
cfg.numrandomization    = 5000;
cfg.design              = design;
cfg.ivar                = 2;
cfg.uvar                = 1;

% run one-sample statistics
stat = ft_freqstatistics(cfg, grand_freq{1}, grand_freq{2});

% extract p-values
[~,~,p] = fdr(squeeze(stat.prob));

% extract cohen's d
for t = 1 : size(stat.stat,3)
    d(t) = computeCohen_d(squeeze(grand_freq{1}.powspctrm(:,1,1,t)),...
                          squeeze(grand_freq{2}.powspctrm(:,1,1,t)),'paired');
end    
     
% calculate number of electrodes showing effect
percent_trending(1,1) = sum((grand_freq{1}.powspctrm(:,1,1,2) - grand_freq{2}.powspctrm(:,1,1,2))<0) ./ size(grand_freq{di}.powspctrm,1);
percent_trending(2,1) = sum((grand_freq{1}.powspctrm(:,1,1,6) - grand_freq{2}.powspctrm(:,1,1,6))>0) ./ size(grand_freq{di}.powspctrm,1);

% save
save([data_dir,'\data\xc\stat_encRetFEContrast.mat'],'p','d','percent_trending')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Run Fixed Effects ANOVA
% load data
load([data_dir,'\data\xc\fixedEffectsPow.mat'])

% predefine vectors for ANOVA analysis
X = []; % power
Y = []; % subject
A = []; % frequency
B = []; % encoding/retrieval

% cycle through both data types
for di = 1 : size(grand_FE,1)
    
    % downsample time domain
    cfg                 = [];
    cfg.win_dur         = 0.1;
    cfg.toi             = [grand_FE{di}.time(1) grand_FE{di}.time(end)];
    tmp1                = sd_downsample_freq(cfg,grand_FE{di,1});
    tmp2                = sd_downsample_freq(cfg,grand_FE{di,2});   
    
    % get hit miss difference
    grand_FE{di,1}.powspctrm = tmp1.powspctrm - tmp2.powspctrm;
    grand_FE{di,1}.time = tmp1.time;    
end

% define encoding and retrieval time index
time_idx = [2 6];

% cycle through each participant
for subj = 1 : size(grand_FE{1,1}.powspctrm)

    % cycle through both data types
    for di = 1 : size(grand_FE,1)

        % extract power (ATL/Gamma/Pre-Post) and add to group vector
        X(end+1:end+2,1) = grand_FE{di,1}.powspctrm(subj,1,[2 3],time_idx(di));
    
        % add frequency for subject
        A(end+1:end+2,1) = [1 2];

        % add trial type for subject
        B(end+1:end+2,1) = [di di];

        % add subject number
        Y(end+1:end+2,1) = [subj subj];
    end    
end

% run 2x2 RM-ANOVA
p = rm_anova2_np(X,Y,A,B,5000);

% combined data into single matrix
E = [X A B];
[~,idx] = sort(A);
E = E(idx,:);

Ei = E(E(:,2) == 1,:);
[~,idx] = sort(Ei(:,3));
Ei = Ei(idx,:);

Ej = E(E(:,2) == 2,:);
[~,idx] = sort(Ej(:,3));
Ej = Ej(idx,:);

E = [Ei;Ej];

% get effect size
e = mes2way(E(:,1),E(:,2:3),'partialeta2');

% save
save([data_dir,'\data\xc\stat_gammaANOVA_FE.mat'],'p','e')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations subj_no file_label


