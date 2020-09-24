function m1_SKTData_psd()


%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1:strfind(codefilepath, 'code') + length('code') - 1);

% add util path
addpath(genpath(fullfile(codefolder, 'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% pipelinefolder
[correspipelinefolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% input setup
folder_input = fullfile(codecorresParentfolder, 'm0_SKTData_extract');

%  animal
[i,j]= regexp(folder_input, 'NHP_\w*');
animal = folder_input(i + length('NHP_'):j);


%% save setup
savefolder = correspipelinefolder;


%% input

% ---parameters for extracting trials from all files
tmin_trial = 1.25;
tmin_reach = 0.5;
tmin_return = 0.5;
t_trial = [-1 1]; % time values respect to target onset


% ---parameters for spectrogram
twin = 0.3;
toverlap = twin * 0.9;
t0  = t_trial(1);

% frequency/time range of intest
f_range_roi = [5 50];
t_roi = t_trial(1);



%% code Start Here

%%%--- all trials extraction ---%%%

% extract all normal trials, nchns * ntemp * ntrials %
pdcond = 'normal';

files = dir(fullfile(folder_input, ['*_' pdcond '*.mat']));
lfptrials_normal = lfptrials_extract_fromfiles(files, t_trial, tmin_trial, tmin_reach, tmin_return);


% extract all mild trials %
pdcond = 'mild';

files = dir(fullfile(folder_input, ['*_' pdcond '*.mat']));
lfptrials_mild = lfptrials_extract_fromfiles(files, t_trial, tmin_trial, tmin_reach, tmin_return);

% extract all moderate trials %
pdcond = 'moderate';

files = dir(fullfile(folder_input, ['*_' pdcond '*.mat']));
lfptrials_moderate = lfptrials_extract_fromfiles(files, t_trial, tmin_trial, tmin_reach, tmin_return);




%%% --- PSD calc and save --- %%%

% extract brainareas
load(fullfile(files(1).folder, files(1).name), 'T_chnsarea');
brainareas = unique(T_chnsarea.brainarea);


for bi = 1: length(brainareas)
    brainarea = brainareas{bi};
    
    % ---- move task
    task = 'move';
    
    if strcmp(task, 'move')
        idx_phase = [500: 1000];
    end
    
    
    if strcmp(task, 'baseline')
        idx_phase = [1: 500];
    end
    
    % extract lfp_normal, lfp_mild, lfp_moderate for psd calc (nchns, ntemp, ntrials)
    eval(['lfp_normal = lfptrials_normal.' brainarea '(:, idx_phase, :);'])
    eval(['lfp_mild = lfptrials_mild.' brainarea '(:, idx_phase, :);'])
    eval(['lfp_moderate = lfptrials_moderate.' brainarea '(:, idx_phase, :);'])
    
    fs = lfptrials_normal.fs;
    
    % plot
    if ~any(strcmp(brainarea, {'STN', 'GP'})) % not DBS
        plotPSD_comp_1chn(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)
    else
        plotPSD_comp_multichns(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)
    end
    clear task idx_phase lfp_* fs
    
    
    % ---- baseline task
    task = 'baseline';
    
    if strcmp(task, 'move')
        idx_phase = [500: 1000];
    end
    
    
    if strcmp(task, 'baseline')
        idx_phase = [1: 500];
    end
    
    % extract lfp_normal, lfp_mild, lfp_moderate for psd calc (nchns, ntemp, ntrials)
    eval(['lfp_normal = lfptrials_normal.' brainarea '(:, idx_phase, :);'])
    eval(['lfp_mild = lfptrials_mild.' brainarea '(:, idx_phase, :);'])
    eval(['lfp_moderate = lfptrials_moderate.' brainarea '(:, idx_phase, :);'])
    
    fs = lfptrials_normal.fs;
    
    % plot
    if ~any(strcmp(brainarea, {'STN', 'GP'})) % not DBS
        plotPSD_comp_1chn(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)
    else
        plotPSD_comp_multichns(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)
    end
    
    clear task idx_phase lfp_* fs
end
close all


function lfp_trials = lfptrials_extract_fromfiles(files, t_trial, tmin_trial, tmin_reach, tmin_return)
%% extract the combined lfp trials from all the files
%
%       input:
%           files = dir(fullfile(inputfolder, ['*_' pdcond '*.mat']));
%
%           t_trial: a vector with two time values respect to target on set e.g[-1 1]
%
%           tmin_trial, tmin_reach, tmin_return: min t_trial, t_reach t_return values;
%                                                ignore the trial less than any of the condition
%
%       return:
%           lfp_trials: extracted all trials struct (e.g lfp_trials.M1: nchns * ntemp * ntrials)


% extract brain areas from files(1)
load(fullfile(files(1).folder, files(1).name), 'T_chnsarea', 'fs');
brainareas = unique(T_chnsarea.brainarea);


% init lfp_trails struct
lfp_trials = struct();
lfp_trials.fs = fs;
nchns_allareas = struct();
for braini = 1: length(brainareas)
    brainarea = brainareas{braini};
    
    eval(['lfp_trials.' brainarea ' = [];']);
    
    eval(['nchns_allareas.' brainarea ' = strcmp(T_chnsarea.brainarea, brainarea);'])
end


nfiles = length(files);
for filei = 1 : nfiles
    
    % load data
    load(fullfile(files(filei).folder, files(filei).name), 'lfpdata', 'fs', 'T_idxevent');
    
    
    nmin_trial = round(tmin_trial * fs);
    nmin_reach = round(tmin_reach * fs);
    nmin_return = round(tmin_return * fs);
    
    
    for triali = 1: height(T_idxevent)
        idx_target = T_idxevent{triali,'TargetTime'};
        idx_touch = T_idxevent{triali, 'TouchTimeix'};
        idx_mouth = T_idxevent{triali, 'MouthTimeix'};
        
        
        if idx_mouth - idx_target > nmin_trial && idx_touch - idx_target > nmin_reach && idx_mouth - idx_touch > nmin_return
            
            n_trial = round(t_trial * fs) + idx_target;
            if n_trial(1) < 1
                n_trial(1) = 1;
            end
            
            % combine to each brain area field
            for braini = 1: length(brainareas)
                brainarea = brainareas{braini};
                
                
                % extract the chns for brainarea
                eval(['chns = nchns_allareas.' brainarea ';'])
                % extract lfp1trial for brainarea
                lfp1trial = lfpdata(chns, n_trial(1): n_trial(2), triali);
                
                
                % combined to particular brainarea field of lfp_trails
                eval(['lfp_trials.' brainarea '  = cat(3, lfp_trials.' brainarea ', lfp1trial);'])
                
                clear brainarea lfp1trial chns
            end
            
            clear ntrial
        end
        
        clear idx_target idx_touch idx_mouth
        
    end
    
    clear lfpdata fs T_idxevent
    clear nmin_trial nmin_reach nmin_return
    clear triali
end


function plotPSD_comp_1chn(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)
%%

% colors setup
color_normal_range = [224, 255, 255] / 255;
color_normal_mean = [0, 0, 255] / 255;
color_mild_range = [255, 228, 225] / 255;
color_mild_mean = [255, 00, 0] / 255;
color_moderate_range = [238, 238, 238] / 255;
color_moderate_mean = [0, 0, 0] / 255;


%%% extract psd_normal %%%
lfp  = lfp_normal;

% averaged across channels and all trials
lfp = squeeze(mean(mean(lfp, 3),1));

%
lfp = zscore(lfp);

% psd using pwelch
[psd, f_psd] = pwelch(lfp, length(lfp), 0, length(lfp), fs);

% normalized psd
psd = (psd - min(psd) )/ (max(psd) - min(psd));

psd_normal = psd;
clear psd lfp


%%% extract psd_mild %%%
lfp  = lfp_mild;

% averaged across channels and all trials
lfp = squeeze(mean(mean(lfp, 3),1));

% psd using pwelch
[psd, ~] = pwelch(lfp, length(lfp), 0, length(lfp), fs);

% normalized psd
psd = (psd - min(psd) )/ (max(psd) - min(psd));

psd_mild = psd;
clear lfp psd



%%% extract psd_moderate %%%
lfp  = lfp_moderate;


% averaged across channels and all trials
lfp = squeeze(mean(mean(lfp, 3),1));

% zscore lfp
lfp = zscore(lfp);

% psd using pwelch
[psd, ~] = pwelch(lfp, length(lfp), 0, length(lfp), fs);

% normalized psd
psd = (psd - min(psd) )/ (max(psd) - min(psd));

psd_moderate = psd;
clear lfp psd




%%% --- plot ---%%%

% index for roi frequency
idx_roi = find(f_psd >= f_range_roi(1) & f_psd <= f_range_roi(2));

% fs_roi, and psd_roi
fs_roi = f_psd(idx_roi);
psd_normal_roi = psd_normal(idx_roi);
psd_mild_roi = psd_mild(idx_roi);
psd_moderate_roi = psd_moderate(idx_roi);


figure;
smooth_span = 5;
plot(fs_roi, smooth(psd_normal_roi, smooth_span), 'DisplayName', 'normal', 'Color', color_normal_mean, 'LineWidth', linewidth)
hold on
plot(fs_roi, smooth(psd_mild_roi,smooth_span), 'DisplayName', 'mild', 'Color', color_mild_mean, 'LineWidth', linewidth)
plot(fs_roi, smooth(psd_moderate_roi, smooth_span), 'DisplayName', 'moderate', 'Color', color_moderate_mean, 'LineWidth', linewidth)

legend()

% title
title([task ' PSD in ' brainarea ])

% save
% save
savefile = fullfile(savefolder, ['psd_' task '_' brainarea]);
saveas(gcf, savefile, 'png');


function plotPSD_comp_multichns(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)
%% lfp_normal, lfp_mild, lfp_moderate: nchns * ntemp * ntrials

% colors setup
color_normal_range = [224, 255, 255] / 255;
color_normal_mean = [0, 0, 255] / 255;
color_mild_range = [255, 228, 225] / 255;
color_mild_mean = [255, 00, 0] / 255;
color_moderate_range = [238, 238, 238] / 255;
color_moderate_mean = [0, 0, 0] / 255;


% averaged all trials lfp: nchns * ntemp
lfpnchns_normal = squeeze(mean(lfp_normal, 3));
lfpnchns_mild = squeeze(mean(lfp_mild, 3));
lfpnchns_moderate = squeeze(mean(lfp_moderate, 3));

nchns = size(lfpnchns_normal, 1);

for chni = 1: nchns
    
    % normal
    lfp = lfpnchns_normal(chni, :);
    
    % zscore lfp
    lfp = zscore(lfp);
    % psd using pwelch
    [psd, f_psd] = pwelch(lfp, length(lfp), 0, length(lfp), fs);
    
    % normalized psd
    psd = (psd - min(psd) )/ (max(psd) - min(psd));
    
    % psd_normal
    psd_normal = psd;
    clear psd lfp
    
    
    
    
    % mild
    lfp = lfpnchns_mild(chni, :);
    
    % zscore lfp
    lfp = zscore(lfp);
    % psd using pwelch
    [psd, ~] = pwelch(lfp, length(lfp), 0, length(lfp), fs);
    
    % normalized psd
    psd = (psd - min(psd) )/ (max(psd) - min(psd));
    
    % psd_normal
    psd_mild = psd;
    clear psd lfp
    
    
    
    % moderate
    lfp = lfpnchns_moderate(chni, :);
    
    % zscore lfp
    lfp = zscore(lfp);
    % psd using pwelch
    [psd, ~] = pwelch(lfp, length(lfp), 0, length(lfp), fs);
    
    % normalized psd
    psd = (psd - min(psd) )/ (max(psd) - min(psd));
    
    % psd_normal
    psd_moderate = psd;
    clear psd lfp
    
    
    
    %%% --- plot ---%%%
    
    % index for roi frequency
    idx_roi = find(f_psd >= f_range_roi(1) & f_psd <= f_range_roi(2));
    
    % fs_roi, and psd_roi
    fs_roi = f_psd(idx_roi);
    psd_normal_roi = psd_normal(idx_roi);
    psd_mild_roi = psd_mild(idx_roi);
    psd_moderate_roi = psd_moderate(idx_roi);
    
    
    figure;
    smooth_span = 5;
    plot(fs_roi, smooth(psd_normal_roi, smooth_span), 'DisplayName', 'normal', 'Color', color_normal_mean, 'LineWidth', linewidth)
    hold on
    plot(fs_roi, smooth(psd_mild_roi,smooth_span), 'DisplayName', 'mild', 'Color', color_mild_mean, 'LineWidth', linewidth)
    plot(fs_roi, smooth(psd_moderate_roi, smooth_span), 'DisplayName', 'moderate', 'Color', color_moderate_mean, 'LineWidth', linewidth)
    
    legend()
    
    % title
    title([task ' PSD in ' brainarea num2str(chni-1) '-' num2str(chni)])
    
    
    % save
    savefile = fullfile(savefolder, ['psd_' task '_' brainarea num2str(chni-1)]);
    saveas(gcf, savefile, 'png');
    
    
    
    clear psd_*  fs_roi f_psd
end





