function m3_SKTData_PlotSpectrogram()
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


%% save setup
savefolder = codecorresfolder;

%% starting: narrow filter the lfp data of all the files
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    eval(['tdur_trial = tdur_trial_' pdcond ';']);
    eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
    eval(['t_minmax_return = t_minmax_return_' pdcond ';']);
    
    
    [lfp_phase_trials, fs, T_chnsarea, idxGroups, idxGroupNames] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach, t_minmax_return);
    plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, idxGroups, idxGroupNames, tdur_trial, fs, animal, pdcond, align2)
    saveas(gcf, fullfile(savefolder, [animal '_' pdcond]), 'png');
    
    clear cond files tdur_trial t_minmax_reach t_minmax_return
    clear lfp_phase_trials fs T_chnsarea
    
    close all
    
end
end


function plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, idxGroups, idxGroupNames, tdur_trial, fs, animal, pdcond, align2)
% plot lfpdata of all the channels: nchns * ntemp * ntrials

twin = 0.2;
toverlap = 0.15;
f_AOI = [5 40];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);


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
nGrps = length(idxGroups);
if (nGrps - 1) * supb_deltaX + nGrps * subp_width +  subp_startLeft > subp_endLeft
    subp_width = round((subp_endLeft - subp_startLeft - (nGrps - 1) * supb_deltaX)/nGrps, 2);
end


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
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    psd_plot = psd(idx_f, :);
    times_plot = times + tdur_trial(1);
    
    % convert into dB and then gauss filted
    psd_plot = 10 * log10(psd_plot);
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[1,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot
end

% plot spectrogram
[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);
ntrials = size(lfp_phase_trials, 3);

fig = figure(); set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
for idxGi = 1 : length(idxGroups)
    idxs = idxGroups{idxGi};
    
    if contains(idxGroupNames{idxGi}, 'STN')
        clim = clim_Spectrogram_STN;
    end
    if contains(idxGroupNames{idxGi}, 'GP')
        clim = clim_Spectrogram_GP;
    end
    if contains(idxGroupNames{idxGi}, 'Others')
        clim = clim_Spectrogram_Others;
    end
    
    
    subp_left = (idxGi -1) * (subp_width + supb_deltaX )+ subp_startLeft;
    
    for idxi = 1 : length(idxs)
        areai = idxs(idxi);
        brainarea = T_chnsarea.brainarea{areai};
        
        subp_bottom = subp_startTop - subp_height - (idxi -1) * (subp_height + supb_deltaY);
        subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
        imagesc(times_plot, freqs_plot, psd_allchns(:, :, idxs(idxi)));
        if isempty(clim)
            set(gca,'YDir','normal')
        else
            set(gca,'YDir','normal', 'CLim', clim)
        end
        
        colormap(jet)
        colorbar
        if idxi == length(idxs)
            xlabel('time/s')
  
            xtls = xticklabels;
            xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
            xticklabels(xtls)
            clear xtls
        else
            xticks([])
        end
        
        ylabel('Frequency(Hz)')
        
        hold on
        % plot reach onset line
        plot([0 0], ylim, 'r--', 'LineWidth',1.5)
        
        
        title([animal ' ' pdcond ':' brainarea ', ntrials = ' num2str(ntrials)])
    end
    
    clear clim
    
end
end
