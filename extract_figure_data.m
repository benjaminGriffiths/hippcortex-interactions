%% Extract Figure Details
% clear workspace
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

%% Figure 1.

%% Figure 2.
% --- resonating freqs --- %
load([data_dir,'\data\res\grand_freq_roi.mat'])
dat = squeeze(grand_roi{1}.powspctrm);
csvwrite([data_dir,'\data\figures\fig2_resHipp.csv'],dat)
dat = squeeze(grand_roi{2}.powspctrm);
csvwrite([data_dir,'\data\figures\fig2_resATL.csv'],dat)
clear dat grand_roi

% --- enc > ret gamma --- %
load([data_dir,'\data\res\grand_freq_gamma.mat'])
dat = squeeze(grand_encoding{1}.powspctrm - grand_retrieval{1}.powspctrm);
csvwrite([data_dir,'\data\figures\fig2_encRetGammaSpecHits.csv'],dat)
dat = squeeze(grand_encoding{2}.powspctrm - grand_retrieval{2}.powspctrm);
csvwrite([data_dir,'\data\figures\fig2_encRetGammaSpecMisses.csv'],dat)
dat = [];
idx = grand_encoding{1}.freq >=40 & grand_encoding{1}.freq <=50;
dat(:,1) = mean(grand_encoding{1}.powspctrm(:,:,idx),3);
dat(:,2) = mean(grand_retrieval{1}.powspctrm(:,:,idx),3);
csvwrite([data_dir,'\data\figures\fig2_45HzGammaRain.csv'],dat)
dat = [];
idx = grand_encoding{1}.freq >=60 & grand_encoding{1}.freq <=80;
dat(:,1) = mean(grand_encoding{1}.powspctrm(:,:,idx),3);
dat(:,2) = mean(grand_retrieval{1}.powspctrm(:,:,idx),3);
csvwrite([data_dir,'\data\figures\fig2_70HzGammaRain.csv'],dat)
clear idx dat grand_encoding grand_retrieval

% --- peak-locked avg. --- %
dat_enc = [];
dat_ret = [];

for subj = 1:nsubj
    
    % load data
    load([data_dir,'\data\preprocessing\pp',num2str(subj),'_data.mat']);
    
    % get peak-locked average
    cfg = [];
    cfg.maxwin = 50;
    cfg.bpfreq = peak_frequencies.enc_gamma(subj,:);
    pf = ft_timelockanalysis([],eeg_getPeakLockedData(cfg,data{1}));
    dat_enc(subj,:) = pf.avg;
    cfg.bpfreq = peak_frequencies.ret_gamma(subj,:);
    pf = ft_timelockanalysis([],eeg_getPeakLockedData(cfg,data{2}));
    dat_ret(subj,:) = pf.avg;        
end
csvwrite([data_dir,'\data\figures\fig2_70HzPeakAvg.csv'],dat_enc)
csvwrite([data_dir,'\data\figures\fig2_45HzPeakAvg.csv'],dat_ret)
clear cfg dat_enc dat_ret pf data subj       

% --- raw trace --- %
load([data_dir,'\data\preprocessing\pp1_data.mat'])
cfg         = [];
cfg.channel = 'Hipp_R_02_4';
cfg.trials  = 75;
cfg.latency = [-1.68 -1.68+0.25];
data        = ft_selectdata(cfg,data{1});
cfg         = [];
cfg.detrend = 'yes';
cfg.hpfilter = 'yes';
cfg.hpfreq   = 1;
cfg.hpfilttype = 'fir';
data        = ft_preprocessing(cfg,data);
dat = data.trial{1};
csvwrite([data_dir,'\data\figures\fig2_70HzRawTrace.csv'],dat)

load([data_dir,'\data\preprocessing\pp1_data.mat'])
cfg         = [];
cfg.channel = 'Hipp_R_02_4';
cfg.trials  = 117;
cfg.latency = [2.17 2.17 + 0.25];
data        = ft_selectdata(cfg,data{2});
cfg         = [];
cfg.detrend = 'yes';
cfg.hpfilter = 'yes';
cfg.hpfreq   = 1;
cfg.hpfilttype = 'fir';
data        = ft_preprocessing(cfg,data);
dat = data.trial{1};
csvwrite([data_dir,'\data\figures\fig2_45HzRawTrace.csv'],dat)
clear cfg data dat

% --- neocortical memory --- %
load([data_dir,'\data\power\grand_encoding.mat'])
dat = squeeze(grand_freq{1}.powspctrm(:,1,2,:) - grand_freq{2}.powspctrm(:,1,2,:));
csvwrite([data_dir,'\data\figures\fig2_atlEncTime.csv'],dat)
dat = [];
idx = grand_freq{1}.time >=0.4 & grand_freq{1}.time <=0.6;
dat(:,1) = mean(grand_freq{1}.powspctrm(:,1,2,idx),4);
dat(:,2) = mean(grand_freq{2}.powspctrm(:,1,2,idx),4);
csvwrite([data_dir,'\data\figures\fig2_atlEncRain.csv'],dat)
load([data_dir,'\data\power\grand_retrieval.mat'])
dat = squeeze(grand_freq{1}.powspctrm(:,1,2,:) - grand_freq{2}.powspctrm(:,1,2,:));
csvwrite([data_dir,'\data\figures\fig2_atlRetTime.csv'],dat)
dat = [];
idx = grand_freq{1}.time >=0.8 & grand_freq{1}.time <=1.2;
dat(:,1) = mean(grand_freq{1}.powspctrm(:,1,2,idx),4);
dat(:,2) = mean(grand_freq{2}.powspctrm(:,1,2,idx),4);
csvwrite([data_dir,'\data\figures\fig2_atlRetRain.csv'],dat)
clear idx dat grand_freq
    
% --- hippocampal ANOVA --- %
dat = [];
for subj = 1 : nsubj

    % load data
    load([data_dir,'\data\power\pp',num2str(subj),'_gamma.mat']);

    % average over channels
    cfg             = [];
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    freq{1}         = ft_selectdata(cfg,freq{1});
    freq{2}         = ft_selectdata(cfg,freq{2});
    
    % add to dat
    dat(subj,1:2) = freq{1}.powspctrm;
    dat(subj,3:4) = freq{2}.powspctrm;
end
csvwrite([data_dir,'\data\figures\fig2_hippGammaANOVA.csv'],dat)
clear dat freq

% --- hippocampal timeseries --- %
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
    clear freq cfg
end

% get grand average
grand_freq          = cell(2,1);
cfg                 = [];
cfg.keepindividual  = 'yes';
for di = 1 : size(group_freq,2)
    grand_freq{di,1} = ft_freqgrandaverage(cfg,group_freq{:,di});
    grand_freq{di,1}.cfg = [];
end
clear cfg group_freq subj di

dat = squeeze(grand_freq{1}.powspctrm(:,:,1,:));
csvwrite([data_dir,'\data\figures\fig2_hippEnc45HzTime.csv'],dat)
dat = squeeze(grand_freq{1}.powspctrm(:,:,2,:));
csvwrite([data_dir,'\data\figures\fig2_hippEnc70HzTime.csv'],dat)
dat = squeeze(grand_freq{2}.powspctrm(:,:,1,:));
csvwrite([data_dir,'\data\figures\fig2_hippRet45HzTime.csv'],dat)
dat = squeeze(grand_freq{2}.powspctrm(:,:,2,:));
csvwrite([data_dir,'\data\figures\fig2_hippRet70HzTime.csv'],dat)
clear dat grand_freq

