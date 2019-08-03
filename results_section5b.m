%% Subsequent Memory Power Analysis
% This script aims to identify the differences in neocortical (NC) and
% medial temporal (MTL) power between later remembered and later forgotten
% pairs.
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
nsubj = numel(dir([data_dir,'\data\preprocessing\pp*']));

%% Get Time-Frequency Data
% cycle through every participant
for subj = 1:nsubj
    
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
    
        % --- get time-frequency data --- %
        % calculate time-frequency
        cfg             = [];
        cfg.channel     = data{di}.label(contact_locations.hippo{subj,1}==1);
        cfg.keeptrials  = 'yes';
        cfg.method      = 'wavelet'; 
        cfg.width       = 5; 
        cfg.output      = 'pow';	
        cfg.pad         = 'nextpow2';
        cfg.foi         = 1.5 : 0.5 : 100;          
        cfg.toi         = toi{di};
        freq{di}       	= ft_freqanalysis(cfg, data{di});

        % subtract 1/f noise
        cfg             = [];
        cfg.toi         = [toi{di}(1) toi{di}(end)];
        freq{di}        = sd_subtr1of(cfg, freq{di});

        % z-transform across frequencies
        freq{di}        = sd_ztransform_freq(freq{di});

        % smooth data
        cfg             = [];
        cfg.fwhm_t      = 0.2;
        cfg.fwhm_f      = 1;
        freq{di}        = smooth_TF_GA(cfg,freq{di});

        % --- split by trial type --- %
        % select hits
        cfg             = [];
        cfg.trials      = freq{di}.trialinfo == 1;
        tmp{1}          = ft_selectdata(cfg,freq{di});
        
        % select misses
        cfg             = [];
        cfg.trials      = freq{di}.trialinfo == 0;
        tmp{2}          = ft_selectdata(cfg,freq{di});
        
        % average over trials
        cfg             = [];
        tmp{1}          = ft_freqdescriptives(cfg, tmp{1});
        tmp{2}          = ft_freqdescriptives(cfg, tmp{2});
        
        % calculate difference
        cfg             = [];
        cfg.parameter   = 'powspctrm';
        cfg.operation   = 'subtract';
        freq{di}        = ft_math(cfg,tmp{1},tmp{2});
                
        % --- split into peak frequencies --- %
        % define frequency fields
        freq_field = {'ret_gamma','enc_gamma'};
        
        % define temporary cell for frequency data
        tmp = cell(numel(freq_field),1);
        
        % cycle through each frequency band
        for fi = 1 : numel(freq_field)

            % select frequencies between sidebands defined in [peak_frequencies] and average
            cfg             = [];
            cfg.frequency   = peak_frequencies.(freq_field{fi})(subj,:);
            cfg.avgoverfreq = 'yes';
            tmp{fi}         = ft_selectdata(cfg, freq{di});
        end

        % append frequencies
        cfg             = [];
        cfg.parameter   = 'powspctrm';
        cfg.appenddim   = 'freq';
        freq{di}        = ft_appendfreq(cfg,tmp{:});

        % re-value frequencies for consistency across participants
        freq{di}.freq   = [45 60];
        freq{di}.time   = linspace(0,1.5,numel(freq{di}.time));
        
        % clear all non-essential variables
        keep home_dir data_dir peak_frequencies contact_locations nsubj subj toi data freq
    end   
    
    % save data
    save([data_dir,'\data\power\pp',num2str(subj),'_gamma.mat'],'freq');
    
    % clear subject data
    clear data toi
end

%% Run Random Effects ANOVA
% predefine vectors for ANOVA analysis
X = []; % power
Y = []; % subject
A = []; % frequency
B = []; % encoding/retrieval

% cycle through each participant
count = 1;
for subj = 1 : nsubj

    % load data
    load([data_dir,'\data\power\pp',num2str(subj),'_gamma.mat']);

    % average over channels
    cfg             = [];
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    freq{1}         = ft_selectdata(cfg,freq{1});
    freq{2}         = ft_selectdata(cfg,freq{2});

    % extract power and add to group vector
    X(end+1:end+2,1) = freq{1}.powspctrm(:);
    X(end+1:end+2,1) = freq{2}.powspctrm(:);
    
    % add frequency for subject
    A(end+1:end+4,1) = [1 2 1 2];
    
    % add trial type for subject
    B(end+1:end+4,1) = [1 1 2 2];
    
    % add subject number
    Y(end+1:end+4,1) = [count count count count];
    count = count + 1;    
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
save([data_dir,'\data\power\stat_gammaRE.mat'],'p','e')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations msubj file_label

%% Run Random Effects T-Tests
% create group cell
group_freq = cell(nsubj,2);

% cycle through each participant
for subj = 1 : nsubj

    % load data
    load([data_dir,'\data\power\pp',num2str(subj),'_gamma.mat']);

    % average over channels
    cfg                 = [];
    cfg.avgoverchan     = 'yes';
    group_freq{subj,1}  = ft_selectdata(cfg,freq{1});
    group_freq{subj,2}  = ft_selectdata(cfg,freq{2});
    
    % relabel channel
    group_freq{subj,1}.label{1} = 'hippo';
    group_freq{subj,2}.label{1} = 'hippo';
end

% get grand average
grand_freq          = cell(2,1);
cfg                 = [];
cfg.keepindividual  = 'yes';
for di = 1 : size(group_freq,2)
    grand_freq{di,1} = ft_freqgrandaverage(cfg,group_freq{:,di});
    grand_freq{di,1}.cfg = [];
end

% define encoding and retrieval file labels
file_label = {'encoding','retrieval'};

% create stat cell
p_fdr	= nan(2,2,7);
d       = nan(2,2,7);
rng(1)

% cycle through each file label
for fi = 1 : numel(file_label)
    
    % downsample time domain
    cfg             = [];
    cfg.win_dur     = 0.2;
    cfg.toi         = [grand_freq{fi}.time(1) grand_freq{fi}.time(end)];
    grand_freq{fi}  = sd_downsample_freq(cfg,grand_freq{fi});
    
    % create null hypothesis
    null_hyp            = grand_freq{fi};
    null_hyp.powspctrm  = zeros(size(null_hyp.powspctrm));
    
    % create design matrix
    design      = [];
    design(1,:) = [1:size(grand_freq{fi}.powspctrm,1), 1:size(grand_freq{fi}.powspctrm,1)];
    design(2,:) = [ones(1,size(grand_freq{fi}.powspctrm,1)) ones(1,size(grand_freq{fi}.powspctrm,1))+1];
    
    figure; hold on
    subplot(1,2,1); hold on 
    plot(grand_freq{fi}.time,squeeze(grand_freq{fi}.powspctrm(:,:,1,:)));
    subplot(1,2,2); hold on 
    plot(grand_freq{fi}.time,squeeze(grand_freq{fi}.powspctrm(:,:,2,:)));
    
    % define default statistic configuration
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'ft_statfun_depsamplesT';
    cfg.correctm            = 'no'; % use fdr function to get adjusted p-values
    cfg.alpha               = 0.05;
    cfg.numrandomization    = 5000;
    cfg.design              = design;
    cfg.ivar                = 2;
    cfg.uvar                = 1;
    
    % --- test hippocampal slow gamma --- %
    cfg.frequency           = [44 46];
    cfg.tail                = 0;
    cfg.channel             = 'hippo';
    stat                    = ft_freqstatistics(cfg,grand_freq{fi},null_hyp);
      
    % extract fdr corrected p-value
    [~,~,p_fdr(fi,1,:)] = fdr(squeeze(stat.prob));
    
    % calculate cohen's d
    for t = 1 : numel(grand_freq{fi}.time)
        d(fi,1,t)           = computeCohen_d(squeeze(grand_freq{fi}.powspctrm(:,1,1,t)),...
                                             squeeze(null_hyp.powspctrm(:,1,1,t)),'paired');
    end
    
    % --- test hippocampal fast gamma --- %
    cfg.frequency           = [59 61];
    cfg.tail                = 0;
    cfg.channel             = 'hippo';
    stat                    = ft_freqstatistics(cfg,grand_freq{fi},null_hyp);
      
    % extract fdr corrected p-value
    [~,~,p_fdr(fi,2,:)] = fdr(squeeze(stat.prob));
    
    % calculate cohen's d
    for t = 1 : numel(grand_freq{fi}.time)
        d(fi,2,t)           = computeCohen_d(squeeze(grand_freq{fi}.powspctrm(:,1,2,t)),...
                                             squeeze(null_hyp.powspctrm(:,1,2,t)),'paired');
    end   
end

% --- read out significance --- %
trl_label  = {'Encoding','Retrieval'};
freq_label = {'45Hz','60Hz'};

% cycle through encoding/retrieval
for ti = 1 : numel(trl_label)
    
    % cycle through channels
    for fi = 1 : numel(freq_label)
        
        % find pfdr < 0.05
        idx = find(p_fdr(ti,fi,:) <= 0.05);
                
        % cycle through each index
        if ~isempty(idx)
            for i = 1:numel(idx)
                fprintf('%s in %s: p = %3.3f, d = %3.3f\n',trl_label{ti},freq_label{fi},p_fdr(ti,fi,idx(i)),abs(d(ti,fi,idx(i))));
            end   
        end
    end
end

% save
save([data_dir,'\data\power\stat_hippocampalGammaRE.mat'],'p_fdr','d')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations nsubj p_fdr d
