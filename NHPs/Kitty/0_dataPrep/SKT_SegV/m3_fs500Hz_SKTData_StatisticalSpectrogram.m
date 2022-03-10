function m3_fs500Hz_SKTData_StatisticalSpectrogram()
%  extract lfp data respect to reachonset
% 
%  return:
%        lfptrials: nchns * ntemp * ntrials


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


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);



%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_goodReach');

cond_cell = cond_cell_extract(animal);
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] ...
    = goodSKTTrials_reachReturn_tcritiria(animal);

if strcmpi(animal, 'bug')
    tdur_trial_normal = [-0.6 1];
    tdur_trial_mild = [-0.6 1];
    tdur_trial_moderate = [-0.6 1];
end
if strcmpi(animal, 'jo') 
    
    tdur_trial_normal = [-0.8 0.8];
    tdur_trial_mild = [-0.8 0.8];
    tdur_trial_moderate = [-0.8 0.8];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    tdur_trial_normal = [-0.6 1];
    tdur_trial_moderate = [-0.6 1];
    
end

if strcmpi(animal, 'pinky')
    tdur_trial_normal = [-0.6 1];
    tdur_trial_mild = [-0.6 1];
    tdur_trial_moderate = [-.6 1];
end

% align to event
align2 = SKTEvent.ReachOnset;


f_AOI = [8 40];
t_AOI = [-0.5 0.5];
twin = 0.2;
toverlap = 0.18;
lowbeta = [8 20];
highbeta = [21 35];
t_preM = [-0.2 0];
t_aftM = [0 0.2];


% saved fig format
savefig_format = 'tif';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

%% Start Code Here

chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    eval(['tdur_trial = tdur_trial_' pdcond ';']);
    eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
    
    
    %%% plot spectrogram across all trials  
    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, t_minmax_reach);
    mask_usedChns = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    lfptrials = lfptrials(mask_usedChns, :, :);
    
    [psd_allchns_alltrials, freqs_used, times_used] = calc_spectrogram_perTrial(lfptrials, tdur_trial, fs, ...
                                  'f_AOI', f_AOI, 't_AOI', t_AOI, 'twin', twin, 'toverlap', toverlap);
                
       
    eval(['psds_' pdcond  '= psd_allchns_alltrials;']);


    clear cond files tdur_trial t_minmax_reach t_minmax_return
    clear lfptrials fs  mask_usedChns
    clear psd_allchns_alltrials
   
end

mask_t_preM = (times_used >= t_preM(1) &times_used <= t_preM(2));
mask_t_aftM = (times_used >= t_aftM(1) &times_used <= t_aftM(2));

betatype = 'lowbeta'
eval(['betabands = ' betatype ';'])
mask_f = (freqs_used >= betabands(1) &freqs_used <= betabands(2));
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    eval(['psd_preM_' pdcond ' = psds_' pdcond '(mask_f, mask_t_preM, :, :);'])
    eval(['psd_aftM_' pdcond ' = psds_' pdcond '(mask_f, mask_t_aftM, :, :);'])    
    clear pdcond
end
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    eval(['avgpsd_preM = squeeze(mean(mean(psd_preM_' pdcond ', 2),1));'])
    eval(['avgpsd_aftM = squeeze(mean(mean(psd_aftM_' pdcond ', 2),1));'])

    for chi = 1 : size(avgpsd_preM,2)
        p = signrank(avgpsd_preM(:, chi),avgpsd_aftM(:, chi));
        
        if p < 0.05
            disp([pdcond ' ' T_chnsarea.brainarea{chi} ', sig diff p = ' num2str(p)])
        end
    end
    clear avgpsd_preM avgpsd_aftM
end
clear betabands mask_f


betatype = 'highbeta'
eval(['betabands = ' betatype ';'])
mask_f = (freqs_used >= betabands(1) &freqs_used <= betabands(2));
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    eval(['psd_preM_' pdcond ' = psds_' pdcond '(mask_f, mask_t_preM, :, :);'])
    eval(['psd_aftM_' pdcond ' = psds_' pdcond '(mask_f, mask_t_aftM, :, :);'])    
    clear pdcond
end
clear betabands mask_f


% statistically 
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    eval(['avgpsd_preM = squeeze(mean(mean(psd_preM_' pdcond ', 2),1));'])
    eval(['avgpsd_aftM = squeeze(mean(mean(psd_aftM_' pdcond ', 2),1));'])

    for chi = 1 : size(avgpsd_preM,2)
        p = signrank(avgpsd_preM(:, chi),avgpsd_aftM(:, chi));
        
        if p < 0.05
            disp([pdcond ' ' T_chnsarea.brainarea{chi} ', sig diff p = ' num2str(p)])
        end
    end
    clear avgpsd_preM avgpsd_aftM
end


end



function [psd_allchns_alltrials, freqs_used, times_used] = calc_spectrogram_perTrial(lfp_phase_trials, tdur_trial, fs, varargin)
% calc spectrogram of lfp_phase_trials: nchns * ntemp * ntrials
%   Inputs:
%       lfp_phase_trials: nchns * ntemp * ntrials
%
%   Returns:
%       psd_allchns: nf * nt  * ntrials * nchns
%       freqs_used: nf * 1
%       times_used: nt * 1


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [-0.5 0.5], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'twin', 0.2, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'toverlap', 0.18, @(x) assert(isnumeric(x) && isscalar(x)));


% parse parameters
parse(p,varargin{:});
f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
twin = p.Results.twin;
toverlap = p.Results.toverlap;


%% Code Start here
nwin = round(twin * fs);
noverlap = round(toverlap * fs);


% calculate psd for each chn across trials
psd_allchns_alltrials = [];
for chi = 1 : size(lfp_phase_trials, 1)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfp_phase_trials, 3)
        x = lfp_phase_trials(chi, :, tri);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_used =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_used = times(idx_t);
    
    psd_plot = psds(idx_f, idx_t, :);
    
    
    psd_allchns_alltrials = cat(4, psd_allchns_alltrials, psd_plot); % psd_allchns: nf * nt * ntrials * nchns
    clear psds  psd_plot idx_t idx_f
end

end
