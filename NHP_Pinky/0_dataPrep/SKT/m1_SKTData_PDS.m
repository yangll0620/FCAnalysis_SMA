function m1_SKTData_PDS()


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
f_range_roi = [10 50];
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

% parameters
brainarea = 'M1';



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
plotPSD_comp_1chn(lfp_normal, lfp_mild, lfp_moderate, F_roi, fs)



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
plotPSD_comp_1chn(lfp_normal, lfp_mild, lfp_moderate, f_range_roi, fs, savefolder, brainarea, task)

