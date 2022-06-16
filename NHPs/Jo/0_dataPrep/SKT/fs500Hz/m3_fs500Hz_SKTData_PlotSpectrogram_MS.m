function m3_fs500Hz_SKTData_PlotSpectrogram_MS()
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


f_AOI = [8 40];
t_AOI = [-0.5 0.5];
twin = 0.2;
toverlap = 0.18;


% saved fig format
savefig_format = 'tif';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

%% Start Code Here

chnsOfI = chnsOfInterest_extract(animal);
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    eval(['tdur_trial = tdur_trial_' pdcond ';']);
    eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
    eval(['t_minmax_return = t_minmax_return_' pdcond ';']);
    
    
    %%% plot spectrogram across all trials  
    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach);
    mask_ChnsOfI = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_ChnsOfI, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    lfptrials = lfptrials(mask_ChnsOfI, :, :);
    
    plot_spectrogram_acrossTrials(lfptrials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, ...
                                  'f_AOI', f_AOI, 't_AOI', t_AOI, 'twin', twin, 'toverlap', toverlap, ... 
                                  'savefolder', savefolder, 'savefig_format', savefig_format);
    saveas(gcf, fullfile(savefolder, [animal '_' pdcond]), savefig_format);
    close gcf
   
    

    %%% final clear
    clear cond files tdur_trial t_minmax_reach t_minmax_return
    clear lfptrials fs T_chnsarea
   
end
end


function plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, varargin)
% plot lfpdata of all the channels: nchns * ntemp * ntrials


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [-0.5 0.5], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'twin', 0.2, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'toverlap', 0.18, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'savefolder', pwd, @isstr);
addParameter(p, 'savefig_format', '.tif', @isstr);


% parse parameters
parse(p,varargin{:});
f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
twin = p.Results.twin;
toverlap = p.Results.toverlap;
savefolder = p.Results.savefolder;
savefig_format = p.Results.savefig_format;



% subplot/Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 1800;
fig_height = 900;


subp_startLeft = 0.05;
subp_endLeft = 0.97;
subp_startTop = 0.98;
subp_width = 0.25;
subp_height = 0.12;
supb_deltaX = 0.02;
supb_deltaY = 0.015;

[psd_allchns, freqs_plot, times_plot] = calc_spectrogram_acrossTrials(lfp_phase_trials, tdur_trial, fs, ...
                                        'f_AOI', f_AOI, 't_AOI', t_AOI, 'twin', twin, 'toverlap', toverlap);

% plot spectrogram
% Group chns into STN, GP and others
mask_STN = contains(T_chnsarea.brainarea, 'stn');
mask_GP = contains(T_chnsarea.brainarea, 'gp');
mask_Others = ~(mask_STN | mask_GP);
idxGroups = [{find(mask_STN)}; {find(mask_GP)}; {find(mask_Others)}];
idxGroupNames = {'STN'; 'GP'; 'Others'};
[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal, 'codesavefolder', fullfile(savefolder, 'code'));
ntrials = size(lfp_phase_trials, 3);

w_map = 200;
h_map = 150;
w_freqNumLabel = 50;
w_colorbar = 50;
h_timeNumLabel = 50;

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
        show_freNumLabel = false;
        show_freqNumLabel = false;
        show_colorbar = false;
        
        fig_width = w_map;
        fig_height = h_map;
        margin_left = 0;
        margin_right = 0;
        margin_top = 0;
        margin_bottom = 0;
        if strcmp(pdcond, 'normal')
            show_freNumLabel = true;
            fig_width = fig_width + w_freqNumLabel;
            margin_left = margin_left + w_freqNumLabel;
        end
        if strcmp(pdcond, 'moderate')
            show_colorbar = true;
            fig_width = fig_width + w_colorbar;
            margin_right = margin_right + w_colorbar;
        end
        if contains(brainarea, 'gp')
            show_freqNumLabel = true;
            fig_height = fig_height + h_timeNumLabel;
            margin_bottom = margin_bottom + h_timeNumLabel;
        end
        
        fig = figure(); 
        set(fig, 'PaperUnits', 'pixels',  'Position', [680 558 fig_width fig_height]);
        ax  = axes('Parent',fig, 'Position', [margin_left margin_bottom fig_width-margin_left-margin_right fig_height-margin_top-margin_bottom]);
        imagesc(ax, times_plot, freqs_plot, psd_allchns(:, :, idxs(idxi)));
        colormap(jet)
        
        if show_freNumLabel
            ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
        else
            yticks([]);
        end
        if show_freqNumLabel
            xlabel('time/s', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
            xtls = xticklabels(ax);
            xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
            xticklabels(ax, xtls)
        else
            xticks([])
        end
        
        
        if show_colorbar
            pos = get(ax, 'Position');
            colorbar('FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
            set(ax, 'Position', pos);
        end

        hold on
        % plot reach onset line
        plot(ax, [0 0], ylim, 'r--', 'LineWidth',1.5)
        if isempty(clim)
            set(ax,'YDir','normal')
        else
            set(ax,'YDir','normal', 'CLim', clim)
        end
        set(ax, 'Position', [0.09 0.15 0.8 0.79])
        
        savefile_sep = fullfile(savefolder, [animal '_sep_' pdcond '_' brainarea]);
        saveas(fig, savefile_sep, savefig_format);
        close(fig)
        
        clear fig_sep ax_sep  savefile_sep  xtls
        
        
        
        clear areai brainarea
    end
    
    clear clim
    
end
end


function [psd_allchns, freqs_plot, times_plot] = calc_spectrogram_acrossTrials(lfp_phase_trials, tdur_trial, fs, varargin)
% calc spectrogram of lfp_phase_trials: nchns * ntemp * ntrials
%   Inputs:
%       lfp_phase_trials: nchns * ntemp * ntrials
%
%   Returns:
%       psd_allchns: nf * nt * nchns
%       freqs_plot: nf * 1
%       times_plot: nt * 1


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

end

