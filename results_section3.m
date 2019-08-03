%% Resonating Frequency Analysis
% This script aims to identify the differences in neocortical (NC) and
% hippocampal resonating frequencies.
clearvars
close all
clc

% define key directory
home_dir = 'E:\bjg335\projects\hippcortex-interactions';
data_dir = 'Y:\projects\intracranial_sync_desync';

% load contact labels
load([home_dir,'\contact_locations.mat'])

% add subfunctions
addpath([home_dir,'\additional_functions'])

% get number of subjects
nsubj = numel(dir([data_dir,'\data\preprocessing\pp*']));

%% Get Time-Frequency Spectrum
% predefine matrix to hold slope values of TF data
beta = [];

% cycle through every participant
for subj = 1 : nsubj
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);

    % pre-define cell for time-frequency data
    freq = cell(size(data));

    % predefine time window of interest based on encoding/retrieval data
    if subj < 8; toi{1} = 3 : 0.025 : 4.5;
    else; toi{1} = 2 : 0.025 : 3.5; end
    toi{2} = 0 : 0.025 : 1.5;
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(data)

        % calculate time-frequency (only on hippocampal channels and recalled trials)
        cfg            = [];
        cfg.channel    = data{di}.label(contact_locations.hippo{subj,1}==1);
        cfg.keeptrials = 'yes';
        cfg.method     = 'wavelet';
        cfg.width      = 5;
        cfg.output     = 'pow';	
        cfg.pad        = 'nextpow2';
        cfg.foi        = 30 : 0.5 : 100;              
        cfg.toi        = toi{di};        
        freq{di,1}     = ft_freqanalysis(cfg, data{di});

        % subtract 1/f noise
        cfg         = [];
        cfg.toi     = [toi{di}(1) toi{di}(end)];
        freq{di,1}  = sd_subtr1of(cfg, freq{di,1});
        
        % extract slope beta
        beta(subj,di,1,:) = mean(mean(freq{di,1}.slope(freq{di,1}.trialinfo == 1,:,:),2),1);
        beta(subj,di,2,:) = mean(mean(freq{di,1}.slope(freq{di,1}.trialinfo >= 0 & freq{di,1}.trialinfo < 1,:,:),2),1);
        
        % smooth data
        cfg         = [];
        cfg.fwhm_t  = 0.2;
        cfg.fwhm_f  = 5;
        freq{di,1}  = smooth_TF_GA(cfg,freq{di,1});

        % rename time bins for consistency between encoding and retrieval
        freq{di,1}.time = linspace(0,1,numel(freq{di,1}.time));
        freq{di,1}.freq = linspace(round(freq{di}.freq(1)),round(freq{di}.freq(end)),numel(freq{di,1}.freq));

        % split data into hits and misses and average
        cfg             = [];
        cfg.trials      = freq{di,1}.trialinfo == 0;
        cfg.avgoverrpt  = 'yes';
        freq{di,2}      = ft_selectdata(cfg,freq{di,1});
        
        cfg             = [];
        cfg.trials      = freq{di,1}.trialinfo == 1;
        cfg.avgoverrpt  = 'yes';
        freq{di,1}      = ft_selectdata(cfg,freq{di,1});
        
        % record slope
        if di == 1; subj_slope = []; end
        subj_slope(:,:,di) = mean(freq{di,1}.powspctrm,3);

        % average over channels and time
        cfg             = [];
        cfg.avgovertime	= 'yes';
        cfg.avgoverchan = 'yes';
        freq{di,1}      = ft_selectdata(cfg,freq{di,1});
        freq{di,2}      = ft_selectdata(cfg,freq{di,2});

        % relabel labels for consistency across participants
        freq{di,1}.label{1} = 'hippo';
        freq{di,2}.label{1} = 'hippo';
        freq{di,1}.cfg      = [];
        freq{di,2}.cfg      = [];
    end
    
    % save data
    save([data_dir,'\data\res\pp',num2str(subj),'_freq_gamma.mat'],'freq');
    fprintf('sub-%02.0f complete...\n',subj)
    
    % plot channelwise difference
    figure('name',sprintf('sub-%02.0f',subj)); hold on
    diff = subj_slope(:,:,1) - subj_slope(:,:,2);
    plot(freq{1}.freq,diff)
    legend(data{1}.label(contact_locations.hippo{subj} == 1));
    
    % clear all non-essential variables
    keep home_dir data_dir contact_locations beta nsubj
end

% save slope data
save([data_dir,'\data\res\slope_betas.mat'],'beta');

%% Get Grand Average
% pre-define cell to hold group data
group_freq = cell(nsubj,2,2);

% cycle through every participant
for subj = 1 : nsubj
    
    % load data
    load([data_dir,'\data\res\pp',num2str(subj),'_freq_gamma.mat'])
           
    % cycle through encoding and retrieval
    for di = 1 : size(freq,2)
        
        % add to group structure
        group_freq{subj,di,1} = freq{di,1}; % add hit data
        group_freq{subj,di,2} = freq{di,2}; % add miss data       
    end
end

% pre-define cell to hold grand-averaged data
grand_encoding  = cell(size(group_freq,3),1);
grand_retrieval = cell(size(group_freq,3),1);

% cycle through hits and misses
for di = 1 : size(group_freq,3)

    % calculate grand average
    cfg                         = [];
    cfg.keepindividual          = 'yes';
    grand_encoding{di,1}        = ft_freqgrandaverage(cfg,group_freq{:,1,di});
    grand_encoding{di,1}.cfg    = [];
    grand_retrieval{di,1}       = ft_freqgrandaverage(cfg,group_freq{:,2,di});
    grand_retrieval{di,1}.cfg   = [];    
end

% save data
save([data_dir,'\data\res\grand_freq_gamma.mat'],'grand_encoding','grand_retrieval');
    
% clear all non-essential variables
keep home_dir data_dir contact_locations nsubj  

%% Run Inferential Statistics
% load data
load([data_dir,'\data\res\grand_freq_gamma.mat'])

% predefine matrices to house p-values and Cohen's d
p_fdr	= nan(2,7);
d       = nan(2,7);
rng(1);

% cycle through hits and misses
for di = 1 : numel(grand_encoding)

    % convert encoding frequency spectrum to time series
    grand_encoding{di,1}.powspctrm    = permute(grand_encoding{di,1}.powspctrm,[1 2 4 3]);
    grand_encoding{di,1}.time         = grand_encoding{di,1}.freq;
    grand_encoding{di,1}.freq         = 2;

    % convert retrieval frequency spectrum to time series
    grand_retrieval{di,1}.powspctrm    = permute(grand_retrieval{di,1}.powspctrm,[1 2 4 3]);
    grand_retrieval{di,1}.time         = grand_retrieval{di,1}.freq;
    grand_retrieval{di,1}.freq         = 2;
    
    % get difference
    grand_diff = grand_encoding{di,1}.powspctrm - grand_retrieval{di,1}.powspctrm;
    shadedErrorBar(grand_encoding{di,1}.time,mean(grand_diff),sem(grand_diff))
    
    % downsample time domain
    cfg                     = [];
    cfg.win_dur             = 10;
    cfg.toi                 = [grand_encoding{di,1}.time(1) grand_encoding{di,1}.time(end)];
    grand_encoding{di,1}    = sd_downsample_freq(cfg,grand_encoding{di,1});
    grand_retrieval{di,1}   = sd_downsample_freq(cfg,grand_retrieval{di,1});           
    
    % create design matrix
    design      = [];
    design(1,:) = [1:size(grand_encoding{di,1}.powspctrm,1), 1:size(grand_encoding{di,1}.powspctrm,1)];
    design(2,:) = [ones(1,size(grand_encoding{di,1}.powspctrm,1)) ones(1,size(grand_encoding{di,1}.powspctrm,1))+1];

    % define statistic configuration
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'ft_statfun_depsamplesT';
    cfg.correctm            = 'no'; % will use fdr function to get corrected p-values
    cfg.correcttail         = 'alpha';
    cfg.alpha               = 0.05;
    cfg.numrandomization    = 5000;
    cfg.design              = design;
    cfg.ivar                = 2;
    cfg.uvar                = 1;

    % run statistics
    stat              = ft_freqstatistics(cfg, grand_encoding{di,1}, grand_retrieval{di,1});
    
    % fdr correct p-values
    [~,~,p_fdr(di,:)]       = fdr(squeeze(stat.prob));
    
    % calculate cohens d
    for t = 1 : numel(grand_encoding{di,1}.time)
        d(di,t) = computeCohen_d(squeeze(grand_encoding{di,1}.powspctrm(:,:,:,t)),squeeze(grand_retrieval{di,1}.powspctrm(:,:,:,t)),'paired');
    end
end

%% Get Peaks
% load data
load([data_dir,'\data\res\grand_freq_gamma.mat'])
load('E:\bjg335\projects\hippcortex-interactions\peak_frequencies.mat')

enc_gamma = nan(nsubj,1);
ret_gamma = nan(nsubj,1);
subj_with_hippo = 1 : nsubj;

for subj = 1:nsubj
        
    figure; hold on

    % get spectrum
    spec = grand_encoding{1}.powspctrm(subj,:) - grand_retrieval{1}.powspctrm(subj,:);
    freq = grand_encoding{1}.freq;
    spec = spec(freq>=50 & freq <= 95);
    freq = freq(freq>=50 & freq <= 95);
    
    % iterate
    success = false;
    itcount = 1;
    fspec = spec;
    while ~success
        
        % find peaks
        [val,idx] = findpeaks(fspec);
        
        % if not empty
        if ~isempty(idx)
            [~,midx] = max(val);
            idx = idx(midx);
            enc_gamma(subj_with_hippo(subj),1) = freq(idx);
            success = true;
            subplot(1,2,1); hold on
            plot(freq,fspec)
            plot(freq(idx),fspec(idx),'o');
            title(num2str(itcount));
        else
            p = polyfit(freq,spec,itcount);
            f = polyval(p,freq);
            fspec = spec - f;
            itcount = itcount + 1;
        end
        if itcount > 3
            success = true;
        end
    end
       
    % get spectrum
    spec = grand_retrieval{1}.powspctrm(subj,:) - grand_encoding{1}.powspctrm(subj,:);
    freq = grand_encoding{1}.freq;
    spec = spec(freq>=32 & freq <= 50);
    freq = freq(freq>=32 & freq <= 50);
    
    % iterate
    success = false;
    itcount = 1;
    fspec = spec;
    while ~success
        
        % find peaks
        [val,idx] = findpeaks(fspec);
        
        % if not empty
        if ~isempty(idx)
            [~,midx] = max(val);
            idx = idx(midx);
            ret_gamma(subj_with_hippo(subj),1) = freq(idx);
            success = true;
            subplot(1,2,2); hold on
            plot(freq,fspec)
            plot(freq(idx),fspec(idx),'o');
            title(num2str(itcount));
        else
            p = polyfit(freq,spec,itcount);
            f = polyval(p,freq);
            fspec = spec - f;
            itcount = itcount + 1;
        end
        if itcount > 3
            success = true;
        end
    end
end

peak_frequencies.enc_gamma = [enc_gamma-5 enc_gamma+5];
peak_frequencies.ret_gamma = [ret_gamma-5 ret_gamma+5];
save peak_frequencies peak_frequencies

%% Run Statistical Analysis on Beta Slopes
% load beta data
load([data_dir,'\data\res\slope_betas.mat'])

% create data
data_enc = struct('cfg',[],...
                  'freq',[1 2],...
                  'label',{{'dummy'}},...
                  'time',1,...
                  'dimord','subj_chan_freq_time',...
                  'powspctrm',mean(beta(:,1,:,:),4));
data_ret = struct('cfg',[],...
                  'freq',[1 2],...
                  'label',{{'dummy'}},...
                  'time',1,...
                  'dimord','subj_chan_freq_time',...
                  'powspctrm',mean(beta(:,2,:,:),4));

% create design matrix
design      = [];
design(1,:) = [1:size(data_ret.powspctrm,1), 1:size(data_ret.powspctrm,1)];
design(2,:) = [ones(1,size(data_ret.powspctrm,1)) ones(1,size(data_ret.powspctrm,1))+1];

% define statistic configuration
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'no'; % will use fdr function to get corrected p-values
cfg.correcttail         = 'alpha';
cfg.alpha               = 0.05;
cfg.numrandomization    = 5000;
cfg.design              = design;
cfg.ivar                = 2;
cfg.uvar                = 1;

% run statistics
rng(1);
stat              = ft_freqstatistics(cfg, data_enc, data_ret);
