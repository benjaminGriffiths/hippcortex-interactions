%% Resonating Frequency Analysis
% This script aims to identify the differences in neocortical (NC) and
% hippocampal resonating frequencies.

clearvars
close all
clc

% define key directory
home_dir = 'E:\bjg335\projects\sync_desync';
data_dir = 'Y:\projects\intracranial_sync_desync';

% load contact labels
load([home_dir,'\contact_locations.mat'])

% add subfunctions
addpath([home_dir,'\additional_functions'])

% define key parameters
n_subj = 8;

%% Get Time-Frequency Spectrum
% cycle through every participant
for subj = 1 : n_subj
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);
    
    % pre-define cell for time-frequency data
    freq = cell(size(data));
    
    % predefine time window of interest based on encoding/retrieval data
    toi{1} = 3 : 0.025 : 4.5;
    toi{2} = 0 : 0.025 : 1.5;
    
    % cycle through encoding and retrieval data
    for di = 1 : numel(data)
    
        % calculate time-frequency
        cfg            = [];
        cfg.keeptrials = 'yes';
        cfg.method     = 'wavelet';
        cfg.width      = 5;
        cfg.output     = 'pow';	
        cfg.pad        = 'nextpow2';
        cfg.foi        = 1 : 0.5 : 100;          
        cfg.toi        = toi{di};
        freq{di}       = ft_freqanalysis(cfg, data{di});

        % subtract 1/f noise
        cfg         = [];
        cfg.toi     = [toi{di}(1) toi{di}(end)];
        freq{di}    = sd_subtr1of(cfg, freq{di});

        % smooth data
        cfg         = [];
        cfg.fwhm_t  = 0.2;
        cfg.fwhm_f  = 1;
        freq{di}    = smooth_TF_GA(cfg,freq{di});
        
        % rename time bins for consistency between encoding and retrieval
        freq{di}.time = linspace(0,1,numel(freq{di}.time));
        freq{di}.freq = linspace(round(freq{di}.freq(1)),round(freq{di}.freq(end)),numel(freq{di}.freq));
    end
    
    % tidy up
    clear cfg di toi data
       
    % concatenate over trials
    cfg            	= [];
    cfg.parameter  	= 'powspctrm';
    freqtmp      	= ft_appendfreq(cfg,freq{:});
    
    % select hippocampus (averaging over contacts, trials and time)
    cfg             = [];
    cfg.avgovertime	= 'yes';
    cfg.avgoverrpt	= 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = freqtmp.label(contact_locations.hippo{subj,1});
    freq_roi{1,1}   = ft_selectdata(cfg,freqtmp);
    
    % select ATL (averaging over contacts, trials and time)
    cfg             = [];
    cfg.avgovertime	= 'yes';
    cfg.avgoverrpt	= 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = freqtmp.label(contact_locations.atl{subj,1});
    freq_roi{2,1}   = ft_selectdata(cfg,freqtmp);
    
    % select PTPR (averaging over contacts, trials and time)
    cfg             = [];
    cfg.avgovertime	= 'yes';
    cfg.avgoverrpt	= 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = freqtmp.label(contact_locations.nc{subj,1});
    freq_roi{3,1}   = ft_selectdata(cfg,freqtmp);
    
    % save data
    save([data_dir,'\data\res\pp',num2str(subj),'_freq_roi.mat'],'freq_roi');
    
    % clear all non-essential variables
    keep n_subj home_dir data_dir contact_locations
end

%% Get Grand Average
% get roi names
rN = {'hippo','atl','ptpr'};

% pre-define cell to hold group data
group_roi = cell(n_subj,numel(rN),2);

% cycle through every participant
for subj = 1 : n_subj
    
    % load data
    load([data_dir,'\data\res\pp',num2str(subj),'_freq_roi.mat'])
    
    % cycle through each roi
    for ri = 1 : size(freq_roi,1)
                
        % add subject data to group cell
        group_roi{subj,ri} = freq_roi{ri,1};

        % rename contact to conform across subject
        group_roi{subj,ri}.label{1} = rN{ri};
    end
end

% pre-define cell to hold grand-averaged data
grand_roi = cell(size(group_roi,2),1);

% cycle through each roi
for ri = 1 : size(grand_roi,1)

    % calculate grand average
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    grand_roi{ri,1}     = ft_freqgrandaverage(cfg,group_roi{:,ri});

    grand_roi{ri,1}.powspctrm = single(grand_roi{ri,1}.powspctrm);
    grand_roi{ri,1}.freq      = single(grand_roi{ri,1}.freq);
    grand_roi{ri,1}.cfg       = [];
end

% save data
save([data_dir,'\data\res\grand_freq_roi.mat'],'grand_roi');
    
% clear all non-essential variables
keep home_dir para
