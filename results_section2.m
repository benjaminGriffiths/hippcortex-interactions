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
n_subj = numel(dir([data_dir,'\data\preprocessing\pp*']));

%% Get Time-Frequency Spectrum
% cycle through every participant
for subj = 1 : n_subj
    
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
    
        % calculate time-frequency
        cfg            = [];
        cfg.keeptrials = 'yes';
        cfg.method     = 'wavelet';
        cfg.width      = 5;
        cfg.output     = 'pow';	
        cfg.pad        = 'nextpow2';
        cfg.foi        = 1.5 : 0.5 : 100;          
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
        freq{di}.freq = linspace(1.5,100,numel(freq{di}.freq));
    end
    
    % tidy up
    clear cfg di toi data
       
    % concatenate over trials
    cfg            	= [];
    cfg.parameter  	= 'powspctrm';
    cfg.appenddim   = 'rpt';
    freqtmp      	= ft_appendfreq(cfg,freq{:});
    
    % plot channels
    pow = squeeze(mean(mean(freqtmp.powspctrm,4),1));
    figure('name',sprintf('sub-%02.0f',subj)); hold on
    subplot(3,1,1); hold on; plot(freqtmp.freq,pow(contact_locations.hippo{subj,1}==1,:)); set(gca,'xscale','log'); legend(freqtmp.label(contact_locations.hippo{subj,1}==1),'location','northeastoutside')
    subplot(3,1,2); hold on; plot(freqtmp.freq,pow(contact_locations.atl{subj,1}==1,:)); set(gca,'xscale','log'); legend(freqtmp.label(contact_locations.atl{subj,1}==1),'location','northeastoutside')
    subplot(3,1,3); hold on; plot(freqtmp.freq,pow(contact_locations.nc{subj,1}==1,:)); set(gca,'xscale','log'); legend(freqtmp.label(contact_locations.nc{subj,1}==1),'location','northeastoutside')
    
    % select hippocampus (averaging over contacts, trials and time)
    cfg             = [];
    cfg.avgovertime	= 'yes';
    cfg.avgoverrpt	= 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = freqtmp.label(find(contact_locations.hippo{subj,1}));
    freq_roi{1,1}   = ft_selectdata(cfg,freqtmp);
    
    % select ATL (averaging over contacts, trials and time)
    cfg             = [];
    cfg.avgovertime	= 'yes';
    cfg.avgoverrpt	= 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = freqtmp.label(find(contact_locations.atl{subj,1}));
    freq_roi{2,1}   = ft_selectdata(cfg,freqtmp);
    
    % select PTPR (averaging over contacts, trials and time)
    cfg             = [];
    cfg.avgovertime	= 'yes';
    cfg.avgoverrpt	= 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = freqtmp.label(find(contact_locations.nc{subj,1}));
    freq_roi{3,1}   = ft_selectdata(cfg,freqtmp);
    
    % save data
    save([data_dir,'\data\res\pp',num2str(subj),'_freq_roi.mat'],'freq_roi');
    fprintf('sub-%02.0f',subj)
    
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
    
%% Plot
figure(); hold on
colmap = [0.8 0.3 0.3; 0.3 0.3 0.8; 0.3 0.8 0.3];
for i = 1 : 3
    x = grand_roi{i}.freq(:,2:end);
    y = squeeze(grand_roi{i}.powspctrm(:,2:end));
    idx = ~any(isnan(y'));
    shadedErrorBar(x,mean(y(idx,:)),sem(y(idx,:)),{'color',colmap(i,:)},1);
end
ylabel('Power (a.u.)'); xlabel('Frequency (Hz)')
xlim([x(1) x(end)])
set(gca,'box','off','tickdir','out','xscale','log',...
    'xtick',[2:2:10 20:20:100]);

%% Get Peaks
theta = nan(n_subj,1);
alpha = nan(n_subj,1);

for subj = 1 : n_subj
        
        figure; hold on

    % get spectrum
    spec = grand_roi{1}.powspctrm(subj,:);
    freq = grand_roi{1}.freq;
    spec = spec(freq>=1.5 & freq <= 7);
    freq = freq(freq>=1.5 & freq<=7);
    
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
            theta(subj,1) = freq(idx);
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
    spec = grand_roi{2}.powspctrm(subj,:);
    freq = grand_roi{2}.freq;
    spec = spec(freq>=7 & freq <= 15);
    freq = freq(freq>=7 & freq<=15);
    
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
            alpha(subj,1) = freq(idx);
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

peak_frequencies.theta = [theta-0.5 theta+0.5];
peak_frequencies.alpha = [alpha-1 alpha+5];
save peak_frequencies peak_frequencies