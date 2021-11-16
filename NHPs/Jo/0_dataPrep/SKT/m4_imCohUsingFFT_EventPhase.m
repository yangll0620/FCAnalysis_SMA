function m4_imCohUsingFFT_EventPhase()
%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);

%% save setup
savefolder = codecorresfolder;

copyfile2folder(codefilepath, savefolder);

%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');

fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;

EventPhases = {'preMove'; 'Anticipated';'earlyReach'; 'Return';'lateReach'};

image_type = 'tif';


twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];


%% Code start here
cond_cell = cond_cell_extract(animal);

EventPhases = SKT_eventPhases_extract();
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal);

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
notAOI_chns = notInterested_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns





for ei = 1: length(EventPhases)
    event = EventPhases{ei};
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
        if isempty(files)
            disp('files is empty!')
            continue;
        end
        
        eval(['t_minmax_reach = t_minmax_reach_' pdcond ';'])
        eval(['t_minmax_return = t_minmax_return_' pdcond ';'])
        
        [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
        
        [ciCoh_NoAmp, f_selected] = ciCohSKT_FFT_NoAmp(lfptrials, fs, f_AOI);
        
        nchns = size(lfptrials,1);
        phis_allchns = [];
        amps_allchns = [];
        for chni = 1: nchns - 1
            for chnj = chni + 1 : nchns
                [phis, amps, ~] = phaseAmp_SKTPerTrial_FFT(lfptrials, fs, f_AOI);
            end
        end
    end
end



