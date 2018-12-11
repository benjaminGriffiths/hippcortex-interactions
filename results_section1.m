%% Behavioural Analysi
clearvars
close all
clc

% define key directory
home_dir = 'E:\bjg335\projects\sync_desync';
data_dir = 'Y:\projects\intracranial_sync_desync';

% add subfunctions
addpath([home_dir,'\additional_functions'])

%% Prepare Data
% predefine vectors
n_modality = zeros(8,3);
n_hit      = zeros(8,3);

% cycle through each subject
for subj = 1 : 8 
    
    % load data
    load([data_dir,'\data\behaviour\pp',num2str(subj),'_beh'])
    
    % cycle through each tria
    for trl = 1 : size(data,1)
        
        % define hit (i.e. correct and confidence is highly confident')
        if data{trl,4}==1 && data{trl,9}==4
            ishit=1;
            
        % define whether miss (i.e. not correct or 'guessed')
        elseif data{trl,4}~=1 || data{trl,9}==1
            ishit=0;
            
        % else, odd trial
        else
            continue
        end
        
        % get trial type
        if strcmpi(data{trl,3},'BIKE') || strcmpi(data{trl,3},'WATERMILL') || strcmpi(data{trl,3},'FARM') || strcmpi(data{trl,3},'UNDERWATER')
            
            % add trial to group data
            n_modality(subj,1)  = n_modality(subj,1) + 1;
            n_hit(subj,1)       = n_hit(subj,1) + ishit;
            
        elseif strcmpi(data{trl,3},'PIANO') || strcmpi(data{trl,3},'GUITAR') || strcmpi(data{trl,3},'ACCORDIAN') || strcmpi(data{trl,3},'TRUMPET')
            
            % add trial to group data
            n_modality(subj,2)  = n_modality(subj,2) + 1;
            n_hit(subj,2)       = n_hit(subj,2) + ishit;
        end
        
        % get grand total
        n_modality(subj,3)  = n_modality(subj,2) + 1;
        n_hit(subj,3)       = n_hit(subj,2) + ishit;
    end           
end

%% Run Analysis
% get mean performance across modalities
mean_performance = nanmean(n_hit(:,3) ./ n_modality(:,3));

% get performance per subject, per modality
percent_per_modality = n_hit ./ n_modality;

% extract non-NaN values
X = percent_per_modality(~isnan(percent_per_modality(:,1)),1);
Y = percent_per_modality(~isnan(percent_per_modality(:,2)),2);

% get mean auditory and visual performance
Xm = mean(X);
Ym = mean(Y);

% run 'true' ttest
[~,~,~,stat] = ttest2(X,Y);
t = stat.tstat;

% concatenate all data points
Z = cat(1,X,Y);

% predefine matrix for permuted t-values
n_perm = 5000;
perm_t = nan(n_perm,1);

% run permutation tstat
for perm = 1 : n_perm
    
    % get random parition in data
    idx = ismember(1:numel(Z),randperm(numel(Z),numel(X)));
    
    % create permutated data
    Xp = Z(idx);
    Yp = Z(~idx);
    
    % run ttest
    [~,~,~,stat] = ttest2(Xp,Yp);
    
    % extract permuted t
    perm_t(perm,1) = stat.tstat;
end

% sort permuted t-values
sorted_t = sort(perm_t);

% find position of observed t in permuted distribution
t_pos = mean(find(t>=sorted_t,1,'last'),find(t<=sorted_t,1,'first'));

% convert position to p-value
p = t_pos ./ n_perm;

% get cohen's d
d = computeCohen_d(X,Y,'independent');
