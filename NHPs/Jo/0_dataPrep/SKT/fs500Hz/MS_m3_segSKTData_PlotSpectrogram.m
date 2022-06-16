function MS_m3_segSKTData_PlotSpectrogram()
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
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '/', '[A-Za-z]*']);
elseif ispc
    % Code to run on Windows platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '\\', '[A-Za-z]*']);
else
    disp('Platform not supported')
end
animal = codecorresfolder(fi + length('NHPs') + 1:j);


%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');

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
coli_align2 = uint32(align2);

% saved fig format
savefig_format = 'tif';

%% save setup
savefolder = codecorresfolder;

%% starting: narrow filter the lfp data of all the files


chnsOfI = chnsOfInterest_extract(animal);
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    if isempty(files)
        continue
    end
    
    eval(['tdur_trial = tdur_trial_' pdcond ';']);
    eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
    eval(['t_minmax_return = t_minmax_return_' pdcond ';']);
    
    
    %%% plot spectrogram across all trials
    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, t_minmax_reach);
    mask_chnsOfI = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_chnsOfI, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    lfptrials = lfptrials(mask_chnsOfI, :, :);
    
    plot_spectrogram_acrossTrials(lfptrials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, savefolder, savefig_format);
    

    %%% final clear
    clear cond files tdur_trial t_minmax_reach t_minmax_return
    clear lfptrials fs T_chnsarea
   
end
end


function plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, varargin)
% plot lfpdata of all the channels: nchns * ntemp * ntrials


if length(varargin) >= 1
    savefolder = varargin{1};
else
    savefolder = pwd;
end

if length(varargin) >= 2
    savefig_format = varargin{2};
else
    savefig_format = 'tif';
end

twin = 0.2;
toverlap = 0.18;
f_AOI = [8 40];
t_AOI = [-0.5 0.5];

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
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end

% plot spectrogram
% Group chns into STN, GP and others
mask_STN = contains(T_chnsarea.brainarea, 'stn');
mask_GP = contains(T_chnsarea.brainarea, 'gp');
mask_Others = ~(mask_STN | mask_GP);
idxGroups = [{find(mask_STN)}; {find(mask_GP)}; {find(mask_Others)}];
idxGroupNames = {'STN'; 'GP'; 'Others'};
[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);
ntrials = size(lfp_phase_trials, 3);

w_map = 300;
h_map = 150;
for idxGi = 1 : length(idxGroups)
    idxs = idxGroups{idxGi};
    
    if strcmpi(idxGroupNames{idxGi}, 'STN')
        clim = clim_Spectrogram_STN;
    end
    if strcmpi(idxGroupNames{idxGi}, 'GP')
        clim = clim_Spectrogram_GP;
    end
    if strcmpi(idxGroupNames{idxGi}, 'Others')
        clim = clim_Spectrogram_Others;
    end
    
   
    
    for idxi = 1 : length(idxs)
        areai = idxs(idxi);
        brainarea = T_chnsarea.brainarea{areai};
        
        
        
        %%%  plot a sepately figure %%%
        fig_sep = figure();
        
        
        set(fig_sep, 'PaperUnits', 'points',  'Position', [680 558 400 200]);
        ax_sep  = axes('Parent',fig_sep);
        set(ax_sep, 'Units', 'pixels','Position', [45 40 290 150])
        
        imagesc(ax_sep, times_plot, freqs_plot, psd_allchns(:, :, idxs(idxi)));
        colormap(jet)
        

        
        if strcmp(pdcond, 'normal')
            ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roma')
        else
            yticks([])
        end
        
        if contains(brainarea, 'gp')
            xlabel('time (s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roma')
            xticks([-0.5 0 0.5])
            xtls = xticklabels(ax_sep);
            xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
            xticklabels(ax_sep, xtls)
        else
            xticks([]);
        end
        
        if strcmp(pdcond, 'moderate')
            pos = get(gca, 'Position');
            c = colorbar;
            set(gca, 'Position', pos);
            clear pos
        end
        
        hold on
        % plot reach onset line
        plot(ax_sep, [0 0], ylim, 'r--', 'LineWidth',1.5)
        if isempty(clim)
            set(ax_sep,'YDir','normal')
        else
            set(ax_sep,'YDir','normal', 'CLim', clim)
        end
        
        %title(ax_sep, [animal ' ' pdcond ':' brainarea ', ntrials = ' num2str(ntrials)]);
        savefile_sep = fullfile(savefolder, [animal '_sep_' pdcond '_' brainarea]);
        saveas(fig_sep, savefile_sep, savefig_format);
        close(fig_sep)
        
        clear fig_sep ax_sep  savefile_sep  xtls
        
        
        
        clear areai brainarea
    end
    
    clear clim
    
end
end

