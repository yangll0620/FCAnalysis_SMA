function m4_SKTData_PlotSpectrogram_freeze2Reach_NoFreezeReach()
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

inputfolder_trials = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');
inputfolder_freezStruct = fullfile(codecorresParentfolder, 'm3_freezeSKTData_EpisodeExtract');

% align to event
align2 = SKTEvent.ReachOnset;




%% save setup
savefolder = codecorresfolder;

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code Start Here
pdcond = 'moderate';
tdur_trial = [-5.5 1.5];
t_plot = [-5 1];
f_plot = [8 40];


useClim = true;
if(isequal(f_plot, [100 400]))
    useClim = false;
end


files_trials = dir(fullfile(inputfolder_trials, ['*_' pdcond '_*.mat']));
files_freezStruct = dir(fullfile(inputfolder_freezStruct, '*.mat'));

chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);


%%% plot spectrogram across all trials
[lfptrials, fs_lfp, matrials, fs_ma,T_chnsarea] = lfp_goodTrials_noFreeze_align2_K(files_trials, files_freezStruct, align2, tdur_trial, 't_min_reach', 0.5);


% chnsofI
ChnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
lfptrials = lfptrials(ChnsOfI_mask, :, :);
T_chnsarea = T_chnsarea(ChnsOfI_mask, :);
clear ChnsOfI_mask

plot_spectrogramMA_acrossTrials(lfptrials, matrials, fs_ma, T_chnsarea, tdur_trial, fs_lfp, savefolder, 't_plot', t_plot, 'f_plot', f_plot, 'useClim', useClim);


%%% final clear
clear cond files_trials tdur_trial t_minmax_reach t_minmax_return
clear lfptrials fs_lfp T_chnsarea

end


function [lfptrials, fs_lfp, matrials, fs_ma,T_chnsarea] = lfp_goodTrials_noFreeze_align2_K(files_trials, files_freezStruct, align2, tdur_trial, varargin)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately
% [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2PeakV(files, [t_AOI(1) t_AOI(2)], 'codesavefolder', savecodefolder);
%
%   not include trials with t_reach <0.2s
%
%         Args:
%             align2: the event to be aligned
%
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%
%
%
%       Name-Value:
%            't_min_reach' - the minimal reach time used, default 0.5 s
%
%
%           't_max_reach' - the max reach time used, default inf
%
%         return:
%             lfptrials: nchns * ntemp * ntrials
%
%             matrials: ntemp * ntrials
%
%             chnAreas:
%
%             fs:


% parse params
p = inputParser;
addParameter(p, 't_min_reach', 0.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 't_max_reach', inf, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
t_min_reach = p.Results.t_min_reach;
t_max_reach = p.Results.t_max_reach;


% code start here
coli_align2 = uint32(align2);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);


load(fullfile(files_trials(1).folder, files_trials(1).name),  'fs_lfp', 'T_chnsarea', 'fs_ma');

freezeFileNames = extractfield(files_freezStruct, 'name');

nfiles = length(files_trials);
lfptrials = [];
matrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files_trials(filei).name;
    load(fullfile(files_trials(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'selectedTrials', 'T_idxevent_ma', 'smoothWspeed_trial');
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    % extract freezeTrials in current file
    bktdt_str = char(regexp(filename, 'bktdt\d+', 'match')); 
    dateofexp_str = char(regexp(filename, '\d{8}', 'match'));
    idx_freeze = find(cellfun(@(x) contains(x, [dateofexp_str '_' bktdt_str]), freezeFileNames));
    befFreezeTrials = []; % trials that have freeze before reach
    if(~isempty(idx_freeze))
        tmp = files_freezStruct(idx_freeze);
        load(fullfile(tmp.folder, tmp.name), 'freezStruct');
        freezEpisodes = freezStruct.freezEpisodes;
        maniFreezeType = freezStruct.optFreezeTypes{end};
        for fri = 1 : length(freezEpisodes)
            if ~strcmp(freezEpisodes{fri}.freezeType, maniFreezeType)
                befFreezeTrials = [befFreezeTrials; freezEpisodes{fri}.triali];
            end
        end
        befFreezeTrials = unique(befFreezeTrials);
        clear tmp freezStruct freezEpisodes
    end
    clear bktdt_str dateofexp_str idx_freeze

    
    ntrials = length(lfpdata);
    for tri = 1: ntrials
        
        % ignore trials marked with 0 and the trial have freeze bef. reach
        if ~selectedTrials(tri) || any(befFreezeTrials == tri)
            continue
        end
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_min_reach  || t_reach > t_max_reach
            clear t_reach
            continue
        end
        
        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
            idx_peakV_ma = idx + idx_reachonset_ma -1;
            t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
            
            if t_reachonset2peakV < 0.3 || t_peakV2reach < 0.3
                clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
                clear t_reachonset2peakV  t_peakV2reach
                continue;
            end
            
            % extract trial with t_dur
            idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
            idx_lfp0 = idx_peakV_lfp;
            idx_ma0 = idx_peakV_ma;
            
            clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
            clear t_reachonset2peakV  t_peakV2reach
            clear idx_peakV_lfp
        else
            idx_lfp0 = T_idxevent_lfp{tri, coli_align2};
            idx_ma0 = T_idxevent_ma{tri, coli_align2};
        end
        
        
        % extract phase for 1 trial
        lfp_1trial = lfpdata{tri};
        idxdur_lfp = round(tdur_trial * fs_lfp) + idx_lfp0;
        idxdur_lfp(1) = idxdur_lfp(1) + 1;
        lfp_phase_1trial = lfp_1trial(:,idxdur_lfp(1) :idxdur_lfp(2));
        
        
        % cat into lfptrials
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        clear idxdur_lfp lfp_1trial lfp_phase_1trial


        % extract ma phase for 1 trial
        ma_1trial = smoothWspeed_trial{tri};
        idxdur_ma = round(tdur_trial * fs_ma) + idx_ma0;
        idxdur_ma(1) = idxdur_ma(1) + 1;
        ma_phase_1trial = ma_1trial(idxdur_ma(1) :idxdur_ma(2), 1);
        
        % cat into lfptrials
        matrials = cat(2, matrials, ma_phase_1trial);
        clear idxdur_lfp ma_1trial ma_phase_1trial

        
        clear t_reach  
    end

    clear befFreezeTrials
    clear filename
    clear('lfpdata', 'T_idxevent_lfp', 'selectedTrials', 'T_idxevent_ma', 'smoothWspeed_trial');
end
end


function plot_spectrogramMA_acrossTrials(lfptrials, matrials, fs_ma, T_chnsarea, tdur_trial, fs_lfp, savefolder, varargin)
% plot lfpdata of all the channels: nchns * ntemp * ntrials
%
% Inputs:
%       lfptrials: nchns * ntemp * ntrials
%
%       matrials: ntemp * ntrials
%   
%       Name-Value:
%            't_plot' - the start and end of plot time, defult [] for all


% parse params
p = inputParser;
addParameter(p, 't_plot', [], @(x) assert(isnumeric(x) && length(x)==2));
addParameter(p, 'f_plot', [], @(x) assert(isnumeric(x) && length(x)==2));
addParameter(p, 'useClim', true, @(x)isscalar(x)&&islogical(x));


parse(p,varargin{:});
t_plot = p.Results.t_plot;
f_plot = p.Results.f_plot;
useClim = p.Results.useClim;


twin = 0.2;
toverlap = 0.18;


nwin = round(twin * fs_lfp);
noverlap = round(toverlap * fs_lfp);

% calculate psd for each chn across trials
psd_allchns = [];
for chi = 1 : size(lfptrials, 1)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfptrials, 3)
        x = lfptrials(chi, :, tri);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs_lfp); % ps: nf * nt
        psds = cat(3, psds, ps);

        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    psd = mean(psds, 3); % psd: nf * nt

    % select freqs/times and corresponding psd
    psd_plot = psd;
    freqs_plot = freqs;
    times_plot = times + tdur_trial(1);
    if ~isempty(f_plot)
        idx_f = (freqs >= f_plot(1) &  freqs <=f_plot(2));
        freqs_plot =  freqs_plot(idx_f);
        psd_plot = psd_plot(idx_f, :);
        clear idx_f
    end

    
    if ~isempty(t_plot)
        idx_t = (times_plot >= t_plot(1) &  times_plot <=t_plot(2));
        times_plot = times_plot(idx_t);
        psd_plot = psd_plot(:, idx_t);
        clear idx_t
    end

    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);

    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end

%
wSpeed = mean(matrials, 2);
times_wSpeed = [1:length(wSpeed)]/fs_ma + tdur_trial(1);


% plot spectrogram
for chi = 1 : height(T_chnsarea)
    brainarea = T_chnsarea.brainarea{chi};
    if contains(brainarea, 'stn')
        brainarea = 'STN';
    end
    if contains(brainarea, 'gp')
        brainarea = 'GP';
    end

    if strcmp(brainarea, 'M1') && useClim
        clim = [-35 -12];
    end
    if  strcmp(brainarea, 'STN') && useClim
        clim = [-30 -6];

    end
    if strcmp(brainarea, 'GP') && useClim
        clim = [-35 -10];
    end

    %%%  plot a sepately figure %%%
    figure(); 

    % plot spectrogram
    subplot(4,1,[2,3,4])
    ax1  = gca;
    imagesc(ax1, times_plot, freqs_plot, psd_allchns(:, :, chi)); hold on;

    if ~exist('clim', 'var')||isempty(clim)
        set(ax1,'YDir','normal')
    else
        set(ax1,'YDir','normal', 'CLim', clim)
    end

    colormap(jet)
    colorbar

    xlabel('time/s')
    ylabel('Frequency(Hz)')
    

    % change time 0 name
    xtls = xticklabels(ax1);
    xtls{find(cellfun(@(x) strcmp(x,'0'), xtls))} = 'ReachOnset';
    xticklabels(ax1, xtls)
    plot(ax1, [0 0], ylim, 'r--', 'LineWidth',1.5)
    pos_ax1 = get(ax1, 'Position');

    % plot wSpeed
    subplot(4,1,1)
    ax2 = gca;
    pos_ax2 = get(ax2, 'Position');
    pos_ax2(3) = pos_ax1(3);
    hlegs = [];
    h = plot(times_wSpeed, wSpeed, 'DisplayName', 'wristSpeed');
    hlegs = [hlegs h]; clear h
    hold on
    h = plot(xlim, [30 30], 'k--', 'DisplayName', 'MoveThres');
    hlegs = [hlegs h];
    clear h
    
    set(ax2, 'Position', pos_ax2);
    xlim(get(ax1, 'XLim'));
    ylimit = ylim;
    if ylimit(2) < 32
        ylim([0 32])
    end    
    plot([0 0], ylim, 'r--', 'LineWidth',1.5)
    set(gca, 'XTick', get(ax1, 'XTick'), 'XTickLabel', get(ax1, 'XTickLabel'))
    ylabel('speed')
    legend(hlegs, 'Position', [0.82 0.85 0.14 0.07])
    clear pos_ax1 pos_ax2 hlegs
    

    pos_ax1 = get(ax1, 'Position');
    pos_ax2 = get(ax2, 'Position');
    pos_ax1(3) = pos_ax2(3);
    set(ax1, 'Position', pos_ax1);

    title(ax2, ['NoFreeze2Reach-MA-spectrogram-' brainarea ' [' num2str(f_plot(1)) ' ' num2str(f_plot(2)) ']Hz'])
    savefile = fullfile(savefolder, ['NoFreeze2Reach-MA-spectrogram-' num2str(f_plot(1)) '-' num2str(f_plot(2)) 'Hz-' brainarea]);
    print(gcf, savefile, '-painters', '-depsc')
    close all
    clear ax1 ax2

end

end

