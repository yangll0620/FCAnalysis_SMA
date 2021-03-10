function m2_SKTData_spectrogram()
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
[fi, j] = regexp(codecorresfolder, fullfile('NHPs','[A-Za-z]*'));
animal = codecorresfolder(fi + length('NHPs/'):j);


%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_avgArea');


t_minmax_reach_normal = [0.7, 1.5];
t_minmax_return_normal = [0.7, 1];
t_minmax_reach_moderate = [0.8, 1.5];
t_minmax_return_moderate = [0.7, 1];

align2 = SKTEvent.ReachOnset;
tdur_trial_normal = [-0.6 1];
tdur_trial_moderate = [-0.6 1];

%% save setup
savefolder = codecorresfolder;

%% starting
files = dir(fullfile(inputfolder, '*.mat'));
for filei = 1: length(files)
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent');
    
    if contains(filename, 'normal')
        tdur_trial = tdur_trial_normal;
        t_minmax_reach = t_minmax_reach_normal;
        t_minmax_return = t_minmax_return_normal;
    else if contains(filename, 'moderate')
            tdur_trial = tdur_trial_moderate;
            t_minmax_reach = t_minmax_reach_moderate;
            t_minmax_return = t_minmax_return_moderate;
        end
    end
    
    ntrials = size(lfpdata, 3);
    for tri = 1: ntrials
        
        % select trials based on reach and return duration
        t_reach = (T_idxevent{tri, coli_reach} - T_idxevent{tri, coli_reachonset}) / fs;
        t_return = (T_idxevent{tri, coli_mouth} - T_idxevent{tri, coli_returnonset}) / fs;
        if t_reach < t_minmax_reach(1) || t_reach > t_minmax_reach(2)
            clear t_reach
            continue
        end
        if t_return < t_minmax_return(1) || t_reach > t_minmax_return(2)
            clear t_return
            continue
        end
        
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs) + T_idxevent{tri, coli_align2};
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
        
    
        
    
    plot_spectrogram(lfptrials, T_chnsarea, fs, animal, cond, align2, tdur_trial)
    
    clear cond files tdur_trial t_minmax_reach t_minmax_return
    clear lfptrials fs T_chnsarea
    
    close all
    
end













function plot_spectrogram(lfptrials, T_chnsarea, fs, animal, pdcond, align2, tdur_trial)
%%
%   Inputs:
%           lfptrials: nchns * ntemp * ntrials
%



% global parameters
twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];


nwin = round(twin * fs);
noverlap = round(toverlap * fs);
nfft = noverlap;

idxSTN = find(strcmp(T_chnsarea.brainarea, 'STN'));
for i = 1 : length(idxSTN)
    T_chnsarea.brainarea{idxSTN(i)} = ['stn' num2str(i-1) '-' num2str(i)];
    
end

idxGP = find(strcmp(T_chnsarea.brainarea, 'GP'));
for i = 1 : length(idxGP)
    T_chnsarea.brainarea{idxGP(i)} = ['gp' num2str(i-1) '-' num2str(i)];
    
end

% plot for each brain area
for areai = 1 : height(T_chnsarea)
    
    brainarea = T_chnsarea.brainarea{areai};
    chnMask =  strcmp(T_chnsarea.brainarea, brainarea);
    
    
    % calculate psds: nf * nt * ntrials
    psds = [];
    for triali = 1: size(lfptrials, 3)
        x = lfptrials(chnMask, : , triali);
        
        [~, f, t, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    psd = mean(psds, 3);
    
    idx_f = (f>=f_AOI(1) &  f<=f_AOI(2));
    f_selected =  f(idx_f);
    psd_selected = psd(idx_f, :);
    t_selected = t + tdur_trial(1);
    
    
    psd_selected = 10 * log10(psd_selected);
    psd_selected = imgaussfilt(psd_selected,'FilterSize',5);
    
    figure
    set(gcf, 'PaperUnits', 'points',  'Position', [675, 550, 700 450]);
    
    % spectrogram subplot
    subplot('Position', [0.1 0.1 0.7 0.8])
    imagesc(t_selected, f_selected, psd_selected);
    colorbar
    set(gca,'YDir','normal')
    
    xlabel('time/s')
    ylabel('psd')
    xticks([-0.5 0 0.5])
    xticklabels({-0.5, char(align2), 0.5})
    
    title([animal ' ' pdcond ': ' brainarea])
    ylim1 = ylim;
    
    
    % psd subplot
    idx_t = (t_selected <0);
    psd_base = mean(psd_selected(:, idx_t), 2);
    psd_phase = mean(psd_selected(:, ~idx_t), 2);
    psd_rel = psd_phase - psd_base;
    subplot('Position', [0.85 0.1 0.1 0.8])
    plot(psd_rel, f_selected);
    ylim(ylim1)
    xlabel('psd diff')
    clear idx_t psd_base psd_phase psd_rel ylim1
    
   
    
    clear brainarea chnMask psds f t idx_f *_selected 
    clear savefile
end






