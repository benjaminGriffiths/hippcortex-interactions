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

% define key parameters
subj_no = 1:8;

% -------------------------------------------------------------- %
% --- All decomposition is reused from scripts for section 4 --- %
% -------------------------------------------------------------- %

%% Run Inferential Statistics
% define encoding and retrieval file labels
file_label = {'encoding','retrieval'};

% create stat cell
p_fdr	= nan(2,7);
d       = nan(2,7);

% cycle through each file label
for fi = 1 : numel(file_label)
    
    % load data
    load([data_dir,'\data\power\grand_',file_label{fi},'.mat'])

    % remove participant 5 (no data)
    grand_freq{1}.powspctrm = grand_freq{1}.powspctrm([1 2 3 4 6 7 8],:,:,:);
    grand_freq{2}.powspctrm = grand_freq{2}.powspctrm([1 2 3 4 6 7 8],:,:,:);
    
    % downsample time domain
    cfg             = [];
    cfg.win_dur     = 0.2;
    cfg.toi         = [grand_freq{1}.time(1) grand_freq{1}.time(end)];
    grand_freq{1}   = sd_downsample_freq(cfg,grand_freq{1});
    grand_freq{2}   = sd_downsample_freq(cfg,grand_freq{2});
    
    % create design matrix
    design      = [];
    design(1,:) = [1:size(grand_freq{1}.powspctrm,1), 1:size(grand_freq{1}.powspctrm,1)];
    design(2,:) = [ones(1,size(grand_freq{1}.powspctrm,1)) ones(1,size(grand_freq{1}.powspctrm,1))+1];
    
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
    
    % --- test atl alpha/beta desync. --- %
    cfg.frequency           = [3 5];
    cfg.tail                = 1;
    cfg.channel             = 'hippo';
    stat                    = ft_freqstatistics(cfg,grand_freq{1},grand_freq{2});
      
    % extract fdr corrected p-value
    [~,~,p_fdr(fi,:)] = fdr(squeeze(stat.prob));
    
    % calculate cohen's d
    for t = 1 : numel(grand_freq{1}.time)
        d(fi,t)           = computeCohen_d(squeeze(grand_freq{1}.powspctrm(:,1,2,t)),...
                                           squeeze(grand_freq{2}.powspctrm(:,1,2,t)),'paired');
    end
end

% --- read out significance --- %
trl_label  = {'Encoding','Retrieval'};

% cycle through encoding/retrieval
for ti = 1 : numel(trl_label)
    
    % find pfdr < 0.05
    idx = find(p_fdr(ti,:) <= 0.05);

    % cycle through each index
    if ~isempty(idx)
        for i = 1:numel(idx)
            fprintf('%s: p = %3.3f, d = %3.3f\n',trl_label{ti},p_fdr(ti,idx(i)),abs(d(ti,idx(i))));
        end   
    end
end

% save
save([data_dir,'\data\power\stat_hippocampalTheta.mat'],'p_fdr','d')

% clear all non-essential variables
keep home_dir data_dir peak_frequencies contact_locations subj_no p_fdr d
