%% Subsequent Memory Power Analysis
% This script aims to identify the differences in neocortical (NC) and
% medial temporal (MTL) power between later remembered and later forgotten
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

% define key parameters (patient 5 has no hippocampal electrodes)
subj_no = [1 2 3 4 6 7 8];

%% Get Time-Frequency Data
% cycle through every participant
for subj = subj_no
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);
    
    % pre-define cell for time-frequency data
    freq = cell(size(data));
    
    % predefine time window of interest based on encoding/retrieval data
    toi{1} = 3 : 0.025 : 4.5;
    toi{2} = 0 : 0.025 : 1.5;
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(data)
    
        % --- get time-frequency data --- %
        % calculate time-frequency
        cfg             = [];
        cfg.channel     = data{di}.label(contact_locations.hippo{subj,1});
        cfg.keeptrials  = 'yes';
        cfg.method      = 'wavelet'; 
        cfg.width       = 5; 
        cfg.output      = 'pow';	
        cfg.pad         = 'nextpow2';
        cfg.foi         = 1 : 0.5 : 100;          
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
        freq_field = {'enc_gamma','ret_gamma'};
        
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
        freq{di}.freq   = [60 45];
        
        % clear all non-essential variables
        keep home_dir data_dir peak_frequencies contact_locations subj_no subj toi data freq
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
for subj = 1 : numel(subj_no)

    % load data
    load([data_dir,'\data\power\pp',num2str(subj_no(subj)),'_gamma.mat']);

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
    Y(end+1:end+4,1) = [subj subj subj subj];
    
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
keep home_dir data_dir peak_frequencies contact_locations subj_no file_label

%% Run Fixed Effects ANOVA
% predefine vectors for ANOVA analysis
X = []; % power
Y = []; % subject
A = []; % frequency
B = []; % encoding/retrieval

% start counter for number of chans
count = 0;

% cycle through each participant
for subj = 1 : numel(subj_no)

    % load data
    load([data_dir,'\data\power\pp',num2str(subj_no(subj)),'_gamma.mat']);

    % average over time
    cfg             = [];
    cfg.avgovertime = 'yes';
    freq{1}         = ft_selectdata(cfg,freq{1});
    freq{2}         = ft_selectdata(cfg,freq{2});

    % get number of channels
    n_chan = numel(freq{1}.label);
    
    % extract power and add to group vector
    X(end+1:end+n_chan*2,1) = freq{1}.powspctrm(:);
    X(end+1:end+n_chan*2,1) = freq{2}.powspctrm(:);
    
    % add frequency for subject
    A(end+1:end+numel(freq{1}.powspctrm)*2,1) = [ones(n_chan,1);ones(n_chan,1)+1;ones(n_chan,1);ones(n_chan,1)+1];
    
    % add trial type for subject
    B(end+1:end+numel(freq{1}.powspctrm)*2,1) = [ones(n_chan*2,1);ones(n_chan*2,1)+1];
    
    % add subject number
    Y(end+1:end+numel(freq{1}.powspctrm)*2,1) = repmat((1:n_chan)'+count,[4 1]);
    count = count + n_chan;
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
save([data_dir,'\data\power\stat_gammaFE.mat'],'p','e')

% --- test population reliability --- %
% cycle through 1000 perms
for perm = 1 : 100
    
    % get resampling indices
    y = datasample(unique(X(:,1)),numel(unique(X(:,1))),'Replace',true);
    
    % create new table
    Xr = [];
    
    % cycle through each resampling index
    for i = 1 : numel(y)
        
        % find parameters for participant y(i)
        Xr(end+1) = X(X(:,1)==y(i),:);
    end
    
    % run 2x2 RM-ANOVA
    [p(perm),f_obs,f_perm] = rm_anova2_np(Xr,Y,A,B,5000);
    
    % acculuate evidence
    f_diff(perm,1) = (f_obs - median(f_perm))>0;
    
    % update command line
    if mod(perm,10)==0
        fprintf('\n%d%% complete...',round(perm/1));
    end
end

% get population reliability
f_reliability = sum(f_diff) ./ numel(f_diff);

% save
save([data_dir,'\data\power\stat_gammaFE_reliability.mat'],'f_reliability')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations subj_no file_label

%% Run Random Effects T-Tests
% create group cell
group_freq = cell(max(subj_no),2);

% cycle through each participant
for subj = subj_no

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
    grand_freq{di,1} = ft_freqgrandaverage(cfg,group_freq{subj_no,di});
    grand_freq{di,1}.cfg = [];
end

% define encoding and retrieval file labels
file_label = {'encoding','retrieval'};

% create stat cell
p_fdr	= nan(2,2,7);
d       = nan(2,2,7);

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
keep home_dir data_dir peak_frequencies contact_locations subj_no p_fdr d

%% Run Random Effects T-Tests
% create group cell
group_freq = {};

% create data template
template = struct('cfg',[],...
                  'freq',[60 45],...
                  'label',{{'hippo'}},...
                  'time',0:0.025:1.5,...
                  'dimord','chan_freq_time');

% cycle through each participant
for subj = subj_no

    % load data
    load([data_dir,'\data\power\pp',num2str(subj),'_gamma.mat']);

    % cycle through each channel
    for chan = 1 : numel(freq{1}.label)
        
        % add template to group freq structure
        group_freq{end+1,1} = template; %#ok<SAGROW>
        group_freq{end,2}   = template;
        
        % add power spectrum to group structure
        group_freq{end,1}.powspctrm = freq{1}.powspctrm(chan,:,:);
        group_freq{end,2}.powspctrm = freq{2}.powspctrm(chan,:,:);       
    end
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
save([data_dir,'\data\power\stat_hippocampalGammaFE.mat'],'p_fdr','d')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations subj_no p_fdr d
