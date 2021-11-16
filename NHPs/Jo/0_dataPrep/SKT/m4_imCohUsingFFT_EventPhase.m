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

ciCohPhasefile_prefix = 'ciCohPhasefile';

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
        
        ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond '_' event '_align2' align2name '.mat']);
        if(~exist(ciCohPhasefile, 'file'))
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            [deltaphis_allChnsTrials, ciCoh, T_chnsarea, ntrials, f_selected]= ciCoh_trialDeltaPhi( inputfolder, pdcond, align2, t_AOI, f_AOI, t_minmax_reach);
            save(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
            clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
            clear t_minmax_reach
        end
        load(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
        
        
        clear pdcond align2 t_AOI align2name ciCohPhasefile
        clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
    end
end



% Histogram chart in polar coordinatescollapse
maxV = max(ciCoh(ciCoh >0));
[chi_maxV, chj_maxV, freqi_maxV] = ind2sub(size(ciCoh), find(ciCoh == maxV));
minV = min(ciCoh(ciCoh >0));
[chi_minV, chj_minV, freqi_minV] = ind2sub(size(ciCoh), find(ciCoh == minV));
deltaPhases_maxV = squeeze(deltaphis_allchns(chi_maxV, chj_maxV, freqi_maxV, :));
deltaPhases_minV = squeeze(deltaphis_allchns(chi_minV, chj_minV, freqi_minV, :));
clear minV chi_minV chj_minV freqi_minV
clear maxV chi_maxV chj_maxV freqi_maxV

nbins = 5;
figure
polarhistogram(deltaPhases_maxV,nbins);
figure
polarhistogram(deltaPhases_minV,nbins);


function [deltaphis_allChnsTrials, ciCoh, T_chnsarea, ntrials, f_selected]= ciCoh_trialDeltaPhi( inputfolder, pdcond, align2, t_AOI, f_AOI, t_minmax_reach)
%
%
%   
%   Output:
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials
%       ciCoh:  nchns * nchns * nf
%       f_selected: nf * 1
%       T_chnsarea: nchns * 5


files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
if isempty(files)
    disp('files is empty!')
    deltaphis_allChnsTrials = [];
    ciCoh = [];
    T_chnsarea = [];
    ntrials = [];
    f_selected = [];
    
    return;
end

[lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);

% extract ciCoh
[ciCoh, f_selected] = ciCohSKT_FFT_NoAmp(lfptrials, fs, f_AOI);

% extract deltaphis for all chns and trials
[nchns, ~, ntrials] = size(lfptrials);
nf = length(f_selected);
deltaphis_allChnsTrials = zeros(nchns, nchns, nf, ntrials);
for chni = 1: nchns - 1
    lfptriali = squeeze(lfptrials(chni, :, :));
    [phisi, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfptriali, fs, f_AOI);
    for chnj = chni + 1 : nchns
        lfptrialj = squeeze(lfptrials(chnj, :, :));
        [phisj, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfptrialj, fs, f_AOI);
        deltaphis_allChnsTrials(chni, chnj, :, :) = phisi - phisj;
        clear lfptrialj phisj
    end
    clear lfptriali phisi
end

