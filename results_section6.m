%% Phase-Amplitude Coupling Analysis
% This script aims to identify the differences in hippocampal theta-gamma 
% phase-amplitude coupling between later remembered and later forgotten
% pairs.
clearvars
close all
clc

% define key directory
home_dir = 'E:\bjg335\projects\sync_desync';
data_dir = 'Y:\projects\intracranial_sync_desync';

% load contact labels
load([home_dir,'\contact_locations.mat'])

% load peak frequencies
load([home_dir,'\peak_frequencies.mat'])

% add subfunctions
addpath([home_dir,'\additional_functions'])

% load experimental parameters
load([data_dir,'\data\para.mat'])

% define key parameters
subj_no     = [1 2 3 4 6 7 8];
time_bins   = {'word','dyn','ret'};
toi         = [3.5 5.5; 0.5 2.5; 0.5 2.5];
time_data   = [1 1 2];

%% Calculate Phase-Amplitude Coupling and Derive Memory Contrast
% predefine cell to contain all PAC data
raw_pac = cell(max(subj_no),size(toi,1));

% cycle through each subject
for subj = subj_no
        
    % pre-define structure for hit > miss phase-amplitude coupling data
    freq = struct('cfg',[],...
                  'label',{{'hippo'}},...
                  'time',1:size(toi,1),...
                  'freq',[45 60],...
                  'powspctrm',[],...
                  'dimord','chan_freq_time');
    
    % cycle through each time window of interest
    for time = 1 : size(toi,1)
        
        % load data
        load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);
        data = data{time_data(time)}; % select data for time window of interest
        
        % select hippocampal channels and restrict time window
        cfg             = [];
        cfg.channel     = data.label(contact_locations.hippo{subj,1});
        cfg.latency     = toi(time,:);
        data            = ft_selectdata(cfg, data);
        
        % pre-define pac matrix
        pac = nan(numel(data.trial),numel(data.label),peak_frequencies.enc_gamma(subj,2)-peak_frequencies.ret_gamma(subj,1)+1,size(peak_frequencies.theta,2));
        
        % cycle through each trial
        for trl = 1 : numel(data.trial)
            
            % extract the hippocampal signal for that trial
            signal         = data.trial{trl};
            
            % define config structure for phase-amplitude coupling analysis
            cfg           = [];
            cfg.Fs        = 500;
            cfg.xlim      = peak_frequencies.theta(subj,:);
            cfg.v         = peak_frequencies.ret_gamma(subj,1) : peak_frequencies.enc_gamma(subj,2);
            cfg.Ncycles   = 5;
            cfg.K         = 1;
            cfg.L         = 0.25;  
            cfg.nFFT      = 1;    
            tmp           = cfc_dist_onetrial(signal,cfg);
            freq_idx      = cfg.v;
                        
            cfg            = [];
            cfg.method     = para.pac.method; 
            tmp            = cfc_quantification(tmp,[],cfg);
            pac(trl,:,:,:) = tmp.CFC;            
        end

        % store all trials for later correlation analysis
        raw_pac{subj,time} = pac;
        
        % get mean phase-amplitude coupling for hits and misses
        avg_hit  = mean(pac(data.trialinfo==1,:,:,:),1);
        avg_miss = mean(pac(data.trialinfo==0,:,:,:),1);
        
        % get difference between hits and misses
        avg_diff = avg_hit - avg_miss;
        
        % collapse over channels and phase-giving frequency
        avg_diff = squeeze(mean(mean(avg_diff,4),2));
        avg_hit  = squeeze(mean(mean(avg_hit,4),2));
            
        % extract indices for encoding- and retrieval-related gamma
        enc_idx = freq_idx >= peak_frequencies.enc_gamma(subj,1) & freq_idx <= peak_frequencies.enc_gamma(subj,2);
        ret_idx = freq_idx >= peak_frequencies.ret_gamma(subj,1) & freq_idx <= peak_frequencies.ret_gamma(subj,2);
        
        % extract memory difference for enc.- and ret.-related gamma
        enc_diff = mean(avg_diff(enc_idx),1);
        ret_diff = mean(avg_diff(ret_idx),1);
        
        % add to freq structure
        freq.powspctrm(:,2,time) = enc_diff;
        freq.powspctrm(:,1,time) = ret_diff;
        
        % add hit data too
        freq.hitspctrm(:,2,time) = mean(avg_hit(enc_idx),1);
        freq.hitspctrm(:,1,time) = mean(avg_hit(ret_idx),1);
        
        % clear non-essential variables
        keep subj subj_no toi time_data para home_dir data_dir peak_frequencies contact_locations freq raw_pac
    end
       
    % save data
    save([data_dir,'\data\pac\pp',num2str(subj),'_data.mat'],'freq');
end

% save data
save([data_dir,'\data\pac\raw_pac.mat'],'raw_pac');
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Calculate Grand Average
group_freq = cell(max(subj_no),1);

% cycle through each subject
for subj = subj_no
    
    % load data
    load([data_dir,'\data\pac\pp',num2str(subj),'_data.mat'])
    
    % add to group cell
    group_freq{subj,1} = freq;
end
    
% get grand average
cfg                 = [];
cfg.keepindividual  = 'yes';
cfg.parameter       = {'powspctrm','hitspctrm'};
grand_freq          = ft_freqgrandaverage(cfg,group_freq{subj_no});

% save data
save([data_dir,'\data\pac\grand_freq.mat'],'grand_freq')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Run Hit > Miss Inferential Statistics
% load data
load([data_dir,'\data\pac\grand_freq.mat'])

% create null hypothesis based on real data
null_hyp            = grand_freq; 
null_hyp.powspctrm  = zeros(size(grand_freq.powspctrm));

% create design matrix
design      = [];
design(1,:) = [1:size(grand_freq.powspctrm,1), 1:size(grand_freq.powspctrm,1)];
design(2,:) = [ones(1,size(grand_freq.powspctrm,1)) ones(1,size(grand_freq.powspctrm,1))+1];

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
cfg.tail                = 1;

% run statistics
stat = ft_freqstatistics(cfg, grand_freq, null_hyp);

% extract p--values
p = squeeze(stat.prob);

% calculate cohen's d
for f = 1 : numel(grand_freq.freq)
    for t = 1 : numel(grand_freq.time)
        d(f,t) = computeCohen_d(squeeze(grand_freq.powspctrm(:,1,f,t)),...
                                squeeze(null_hyp.powspctrm(:,1,f,t)),'paired');
    end   
end

% save
save([data_dir,'\data\pac\stat_hitMissContrast.mat'],'p','d')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Contrast Verbal and Dynamic SMEs
% load data
load([data_dir,'\data\pac\grand_freq.mat'])

% split data
cfg = [];
cfg.latency = [1 1];
verbPAC = ft_selectdata(cfg,grand_freq);

cfg = [];
cfg.latency = [2 2];
dynPAC  = ft_selectdata(cfg,grand_freq);
dynPAC.time = 1;

% create design matrix
design      = [];
design(1,:) = [1:size(grand_freq.powspctrm,1), 1:size(grand_freq.powspctrm,1)];
design(2,:) = [ones(1,size(grand_freq.powspctrm,1)) ones(1,size(grand_freq.powspctrm,1))+1];

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
cfg.tail                = 1;

% run statistics
stat = ft_freqstatistics(cfg, verbPAC, dynPAC);

% extract p--values
p = squeeze(stat.prob);

% calculate cohen's d
for f = 1 : numel(grand_freq.freq)
    d(f,1) = computeCohen_d(squeeze(verbPAC.powspctrm(:,1,f)),...
                            squeeze(dynPAC.powspctrm(:,1,f)),'paired');
end

% save
save([data_dir,'\data\pac\stat_verbDynContrast.mat'],'p','d')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Contrast Verbal Fast and Slow SMEs
% load data
load([data_dir,'\data\pac\grand_freq.mat'])

% split data
cfg = [];
cfg.latency = [1 1];
cfg.frequency = [59 61];
fastPAC = ft_selectdata(cfg,grand_freq);

cfg = [];
cfg.latency = [1 1];
cfg.frequency = [44 46];
slowPAC = ft_selectdata(cfg,grand_freq);
slowPAC.freq = 60;

% create design matrix
design      = [];
design(1,:) = [1:size(grand_freq.powspctrm,1), 1:size(grand_freq.powspctrm,1)];
design(2,:) = [ones(1,size(grand_freq.powspctrm,1)) ones(1,size(grand_freq.powspctrm,1))+1];

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
cfg.tail                = 1;

% run statistics
stat = ft_freqstatistics(cfg, fastPAC, slowPAC);

% extract p--values
p = squeeze(stat.prob);

% calculate cohen's d
d = computeCohen_d(squeeze(fastPAC.powspctrm),...
                   squeeze(slowPAC.powspctrm),'paired');

% save
save([data_dir,'\data\pac\stat_fastSlowContrast.mat'],'p','d')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no

%% Contrast Verbal SME and RSE
% load data
load([data_dir,'\data\pac\grand_freq.mat'])

% split data
cfg = [];
cfg.latency = [1 1];
encPAC = ft_selectdata(cfg,grand_freq);

cfg = [];
cfg.latency = [3 3];
retPAC  = ft_selectdata(cfg,grand_freq);
retPAC.time = 1;

% create design matrix
design      = [];
design(1,:) = [1:size(grand_freq.powspctrm,1), 1:size(grand_freq.powspctrm,1)];
design(2,:) = [ones(1,size(grand_freq.powspctrm,1)) ones(1,size(grand_freq.powspctrm,1))+1];

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
cfg.tail                = 1;

% run statistics
stat = ft_freqstatistics(cfg, encPAC, retPAC);

% extract p--values
p = squeeze(stat.prob);

% calculate cohen's d
for f = 1 : numel(grand_freq.freq)
    d(f,1) = computeCohen_d(squeeze(encPAC.powspctrm(:,1,f)),...
                            squeeze(retPAC.powspctrm(:,1,f)),'paired');
end

% select fast gamma for encoding and slow gamma for retrieval
cfg = [];
cfg.frequency = [59 61];
encPAC = ft_selectdata(cfg,encPAC);

cfg = [];
cfg.frequency = [44 46];
retPAC  = ft_selectdata(cfg,retPAC);
retPAC.freq = 60;

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
cfg.tail                = 1;

% run statistics
stat = ft_freqstatistics(cfg, encPAC, retPAC);

% extract p- and t-values
p(3) = stat.prob;
d(3) = computeCohen_d(squeeze(encPAC.powspctrm),...
                      squeeze(retPAC.powspctrm),'paired');

% save
save([data_dir,'\data\pac\stat_encRetContrast.mat'],'p','d')

% clear all non-essential variables
keep contact_locations data_dir home_dir peak_frequencies subj_no
