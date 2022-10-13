function m3_SKTData_PSD()


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
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');

cond_cell = {'normal', 'PD'};
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] ...
    = goodSKTTrials_reachReturn_tcritiria(animal);


% align to event
align2 = SKTEvent.ReachOnset;


f_AOI = [8 40];
t_AOI = [];
twin = 0.2;
toverlap = 0.18;


% saved fig format
savefig_format = 'tif';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

%% Start Code Here
EventPhases = {'preMove';'earlyReach'};
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
tdur_trial = [-1, 1];

for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    if strcmpi(pdcond, 'PD')
        files1 = dir(fullfile(inputfolder, ['*_mild_*.mat']));
        files2 = dir(fullfile(inputfolder, ['*_moderate_*.mat']));
        files = [files1; files2];
        t_minmax_reach = 0.6;
        clear files1 files2
    else
        files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
        t_minmax_reach = 0.5;
    end
    
    
    %%% plot spectrogram across all trials
    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach);
    mask_ChnsOfI = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_ChnsOfI, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    lfptrials = lfptrials(mask_ChnsOfI, :, :);
    
    [psd_allchns, freqs_plot, times_plot] = calc_spectrogram_acrossTrials(lfptrials, tdur_trial, fs, ...
        'f_AOI', f_AOI, 't_AOI', t_AOI, 'twin', twin, 'toverlap', toverlap);
    
    f_selected = freqs_plot;
    brainareas = {};
    for bi = 1 : height(T_chnsarea)
        brainarea = T_chnsarea.brainarea{bi};
        if(contains(brainarea, 'stn'))
            brainarea = 'STN';
        end
        if(contains(brainarea, 'gp'))
            brainarea = 'GP';
        end
        brainareas = [brainareas; brainarea];
        clear brainarea
    end
    for ei = 1: length(EventPhases)
        eventname = EventPhases{ei};
        if(strcmpi(eventname, 'preMove'))
            idx_t = (times_plot>=-1 & times_plot<=-0.8);
        end
        if(strcmpi(eventname, 'earlyReach'))
            idx_t = (times_plot>=0 & times_plot<=0.2);
        end
        tmp = squeeze(mean(psd_allchns(:, idx_t, :, :), 2)); % tmp: nf * ntrials * nchns
        for bi = 1 : length(brainareas)
            brainarea = brainareas{bi};
            pxx.(brainarea) = squeeze(tmp(:, :, bi));
            clear brainarea
        end
        
        pxxs.(eventname).(pdcond) = pxx;
        clear eventname pxx
    end
    
    %%% final clear
    clear cond files t_minmax_reach t_minmax_return
    clear lfptrials fs T_chnsarea
    
end

save(fullfile(savefolder, [animal '_psd.mat']), 'pxxs', 'f_selected', 'brainareas');
end


function [psd_allchns, freqs_plot, times_plot] = calc_spectrogram_acrossTrials(lfp_phase_trials, tdur_trial, fs, varargin)
% calc spectrogram of lfp_phase_trials: nchns * ntemp * ntrials
%   Inputs:
%       lfp_phase_trials: nchns * ntemp * ntrials
%
%   Returns:
%       psd_allchns: nf * nt * ntrials * nchns
%       freqs_plot: nf * 1
%       times_plot: nt * 1


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [], @(x) assert(isempty(x)||(isnumeric(x) && isvector(x) && length(x)==2)));
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
psd_allchns = [];
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
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    if ~isempty(t_AOI)
        idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
        times_plot = times(idx_t);
        psd_plot = psds(idx_f, idx_t, :);
        clear idx_t
    else
        times_plot = times;
        psd_plot = psds(idx_f, :, :);
    end
    
    
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(4, psd_allchns, psd_plot); % psd_allchns: nf * nt * ntrials * nchns
    clear psds psd psd_plot idx_f
end

end