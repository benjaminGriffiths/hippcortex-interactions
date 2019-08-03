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

% get number of subjects
nsubj = numel(dir([data_dir,'\data\preprocessing\pp*']));

%% Get ERP
% cycle through every participant
for subj = 1:nsubj
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat'],'data');
    
    % pre-define cell for time-frequency data
    freq = cell(size(data));
    
    % predefine time window of interest based on encoding/retrieval data
    if subj < 8; toi{1} = [2.5 4];
    else; toi{1} = [1.5 3]; end
    toi{2} = [-0.5 1];
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(data)
        
        % filter data
        cfg = [];
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 20;
        data{di} = ft_preprocessing(cfg,data{di});
        
        % get timelock
        cfg         = [];
        cfg.trials  = data{di}.trialinfo == 1;        
        cfg.latency = toi{di};
        tml_h       = ft_timelockanalysis(cfg,data{di});
        cfg.trials  = data{di}.trialinfo == 0;    
        tml_m       = ft_timelockanalysis(cfg,data{di});
        
        % baseline
        cfg = [];
        cfg.baseline = [0.25 0.5]+tml_h.time(1);
        tml_h = ft_timelockbaseline(cfg,tml_h);
        tml_m = ft_timelockbaseline(cfg,tml_m);
        
        % get difference
        tml_diff = tml_h;
        tml_diff.avg = tml_h.avg - tml_m.avg;
         
        % --- split by channels --- %
        % define roi fields
        chan_field = fieldnames(contact_locations);
        
        % define temporary cell for frequency data
        tmp = cell(numel(chan_field),1);
        
        % cycle through each roi
        for ci = 1 : numel(chan_field)

            % select frequencies between sidebands define in [para] and average
            cfg             = [];
            cfg.channel     = tml_diff.label(contact_locations.(chan_field{ci}){subj,1}==1);
            cfg.avgoverchan = 'yes';
            tmp{ci,1}      = ft_selectdata(cfg, tml_diff);
        end

        % append channels
        cfg             = [];
        cfg.appenddim   = 'chan';
        tml_diff        = ft_appenddata(cfg,tmp{:});

        % retimelock
        tml_diff        = ft_timelockanalysis([],tml_diff);
        tml_diff.label  = chan_field;
        tml_diff.time   = linspace(-0.5,1,numel(tml_diff.time));
        tml_diff.cfg    = [];      
        
        % add to group
        group_tml{subj,di} = tml_diff;
    end
end

% get grand average
cfg = [];
cfg.keepindividual = 'yes';
grand_encoding = ft_timelockgrandaverage(cfg,group_tml{:,1});
grand_retrieval = ft_timelockgrandaverage(cfg,group_tml{:,2});

% downsample time domain
cfg             = [];
cfg.win_dur     = 0.2;
cfg.toi         = [grand_encoding.time(1) grand_encoding.time(end)];
grand_encoding   = sd_downsample_freq(cfg,grand_encoding);
grand_retrieval  = sd_downsample_freq(cfg,grand_retrieval);

% create null
null_hyp = grand_encoding;
null_hyp.individual = zeros(size(null_hyp.individual));

% create design matrix
design      = [];
design(1,:) = [1:size(grand_encoding.individual,1), 1:size(grand_encoding.individual,1)];
design(2,:) = [ones(1,size(grand_encoding.individual,1)) ones(1,size(grand_encoding.individual,1))+1];

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
stat.enc                = ft_timelockstatistics(cfg,grand_encoding,null_hyp);
stat.ret                = ft_timelockstatistics(cfg,grand_retrieval,null_hyp);

% cycle through each channel
for chan = 1 : 3

    % extract fdr corrected p-value
    [~,~,p_fdr(1,chan,:)] = fdr(squeeze(stat.enc.prob(chan,:)));
    [~,~,p_fdr(2,chan,:)] = fdr(squeeze(stat.ret.prob(chan,:)));
end