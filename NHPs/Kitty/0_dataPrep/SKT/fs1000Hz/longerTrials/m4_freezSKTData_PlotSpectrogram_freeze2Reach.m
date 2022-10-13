function m4_freezSKTData_PlotSpectrogram_freeze2Reach()
%
%
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

pdcond = 'moderate';

%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm3_freezeSKTData_EpisodeExtract');


f_AOI = [8 40];
useClim = true;
if(isequal(f_AOI, [100 150]))
    useClim = false;
end

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

%% Code Start Here
files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
if isempty(files)
    return;
end
    
    
% extract freeze tbl_freezEpisodes
[tbl_freezEpisodes]= tbl_freezEpisodes_extract(files, 't_bef', 1.5, 't_aft', 1.5);
tbl_reachFreezEpisodes = tbl_freezEpisodes(strcmp(tbl_freezEpisodes.freezType, "reachFreeze"),:);
clear tbl_freezEpisodes


plot_spectrogramMA_freeze2Reach_acrossTrials(tbl_reachFreezEpisodes, f_AOI, savefolder, 'useClim', useClim);



end


function [tbl_freezEpisodes]= tbl_freezEpisodes_extract(files, varargin)
% extract tbl_freezEpisodes from all files, not include the bad trials
%
% 
%   Inputs
%       Name-Value: 
%           't_bef': including t_bef before freeze start time (default = 0.8),
%
%           't_aft': including t_aft after freeze end time (default = 1),

if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    tbl_freezEpisodes = [];
    
    return;
end


p = inputParser;
addParameter(p, 't_bef',0.8, @(x)isscalar(x)&&isnumeric(x));
addParameter(p, 't_aft',1, @(x)isscalar(x)&&isnumeric(x));
parse(p,varargin{:});
t_bef = p.Results.t_bef;
t_aft = p.Results.t_aft;


[~, combFreeTypes] = optFreezeTypes_extract();


tblvarNames = {'datebk_str', 'freezi', 'freezType', 'triali', 'lfp', 'idx_lfpStrendFreez', 'fs_lfp', 'T_chnsarea', 'wSpeed', 'idx_wSpeedStrendFreez', 'fs_ma'};


for fi = 1: length(files)
    loadfilename = files(fi).name;
    load(fullfile(files(fi).folder, loadfilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'freezStruct', 'selectedTrials', 'smoothWspeed_trial', 'fs_ma');
    
    datebk_str = regexp(loadfilename, '\d{8}_bktdt\d*', 'match');
    datebk_str = datebk_str{1};
    
    % check fs and T_chnsarea same across files
    if ~exist('fs_unit', 'var')
        fs_unit = fs_lfp;
    else
        if(fs_unit ~= fs_lfp)
            disp(['fs ~= fs_lfp for file ' loadfilename])
            continue;
        end
    end
    if ~exist('T_chnsarea_unit', 'var')
        T_chnsarea_unit = T_chnsarea;
    else
        if(~isequal(T_chnsarea_unit, T_chnsarea))
            disp(['T_chnsarea_unit ~= T_chnsarea for file ' loadfilename])
            continue;
        end
    end
    
    freezEpisodes = freezStruct.freezEpisodes;
    for frzi = 1 : length(freezEpisodes) 
        
        %%%  extract freeze lfp data: seglfp
        tri = freezEpisodes{frzi}.triali;
        if ~selectedTrials(tri)
            clear tri
            continue;
        end
        
        % extract lfp data in freeze phase (including t_bef to t_aft)
        n_bef = round(t_bef * fs_lfp);
        n_aft = round(t_aft * fs_lfp);
        idx_durFreez = round(freezEpisodes{frzi}.freezeTPhaseS * fs_lfp);
        idx_dur = [idx_durFreez(1) - n_bef idx_durFreez(2) + n_aft];
        lfpfreez  = lfpdata{tri}(:, idx_dur(1): idx_dur(2));
        idx_lfpfreezeInExtract = [n_bef+1 length(lfpfreez)-n_aft];
        clear n_bef n_aft idx_durFreez idx_dur
        
        % extract ma data in freeze phase (including t_bef to t_aft)
        n_bef = round(t_bef * fs_ma);
        n_aft = round(t_aft * fs_ma);
        idx_durFreez = round(freezEpisodes{frzi}.freezeTPhaseS * fs_ma);
        idx_dur = [idx_durFreez(1) - n_bef idx_durFreez(2) + n_aft];
        wSpeedfreez  = smoothWspeed_trial{tri}(idx_dur(1): idx_dur(2), :);
        idx_wSpeedfreezeInExtract = [n_bef+1 length(wSpeedfreez)-n_aft];
        clear n_bef n_aft idx_durFreez idx_dur
        
        
        % extract freezType
        idx_freez = find(cellfun(@(x) contains(freezEpisodes{frzi}.freezeType, x), {'init', 'Reach', 'Manipulation'}));  
        freezType = combFreeTypes{idx_freez};
        clear idx_freeze

        
        %%% append to tbl_freezEpisodes
        tbl = table(string(datebk_str), frzi, string(freezType), tri, {lfpfreez}, {idx_lfpfreezeInExtract}, fs_lfp, {T_chnsarea}, {wSpeedfreez}, {idx_wSpeedfreezeInExtract}, fs_ma, 'VariableNames', tblvarNames);
        if ~exist('tbl_freezEpisodes', 'var')
            tbl_freezEpisodes = tbl;
        else
            tbl_freezEpisodes = [tbl_freezEpisodes; tbl];
        end
        clear tbl tri freezType lfpfreez idx_freezeInExtract
       
    end   
    clear freezEpisodes
    clear('lfpdata', 'fs_lfp', 'T_chnsarea', 'freezStruct', 'selectedTrials');
end
end

function plot_spectrogramMA_freeze2Reach_acrossTrials(tbl_freezEpisodes, f_AOI, savefolder, varargin)
% plot spectrogram and MA


p = inputParser;
addParameter(p, 'useClim', true, @(x)isscalar(x)&&islogical(x));
parse(p,varargin{:});
useClim = p.Results.useClim;


if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end



time0name = 'freezeEnd';
tdur_trial = [-5.2 1];
t_plot = [-5 1];


% calc each freeze phase
psds_alltrials_AOI = []; % nf * nt * nchns * ntrials
wSpeeds = []; % wSpeeds: ntemp * ntrials
for tbi = 1 : height(tbl_freezEpisodes)
    
    % extract lfp_freeze duration
    lfp = tbl_freezEpisodes.lfp{tbi}; % lfp: nchns * ntemp
    fs_lfp = tbl_freezEpisodes.fs_lfp(tbi);

    idx_lfpStrendFreez = tbl_freezEpisodes.idx_lfpStrendFreez{tbi};
    
    t_freeze = (idx_lfpStrendFreez(2) - idx_lfpStrendFreez(1))/fs_lfp;
    if t_freeze < abs(tdur_trial(1))
        disp(['less than ' num2str(abs(tdur_trial(1)))])
        continue;
    end
    if strcmp(time0name, 'freezeEnd')
        idx_lfp0 = idx_lfpStrendFreez(2);
    else
        idx_lfp0 = idx_lfpStrendFreez(1);
    end
    idx_lfp = idx_lfp0 + round(fs_lfp * tdur_trial);
    lfp = lfp(:,idx_lfp(1): idx_lfp(2));
    clear idx_lfp
    
    % extract wSpeed_freeze duration
    fs_ma = tbl_freezEpisodes.fs_ma(tbi);
    idx_wSpeedStrendFreez = tbl_freezEpisodes.idx_wSpeedStrendFreez{tbi};
    if strcmp(time0name, 'freezeEnd')
        idx_ma0 = idx_wSpeedStrendFreez(2);
    else
        idx_ma0 = idx_wSpeedStrendFreez(1);
    end
    idx_ma = idx_ma0 + round(fs_ma * tdur_trial);
    wSpeed = tbl_freezEpisodes.wSpeed{tbi};
    wSpeed = wSpeed(idx_ma(1): idx_ma(2), 1); % wSpeed: ntemp * 1
    wSpeeds = cat(2, wSpeeds, wSpeed);
    clear fs_ma  idx_wSpeedStrendFreez wSpeed
    clear idx_ma0 idx_ma
    
    % calc psd
    [psd, freqs, times] = calc_psd_allchns(lfp, fs_lfp);% psds: nf * nt * nchns
    
    if tbi == 1
        % extract based on f_AOI and align times with t_lfpfreezstr as time 0s
        idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        freqs_plot =  freqs(idx_f);
        times_plot = times + tdur_trial(1);
        if ~isempty(t_plot)
            idx_t = (times_plot >= t_plot(1) &  times_plot <=t_plot(2));
            times_plot = times_plot(idx_t);
        end
        clear idx_lfpfreezstr   
    end

    if ~isempty(t_plot)
        psd_AOI = psd(idx_f, idx_t, :);
    else
        psd_AOI = psd(idx_f, :, :);
    end
    
    % append
    psds_alltrials_AOI = cat(4,psds_alltrials_AOI,psd_AOI);
   
    clear psds lfp idx_lfpStrendFreez psd_AOI freqs times fs_lfp
end
psds_AOI = mean(psds_alltrials_AOI, 4);


wSpeed = mean(wSpeeds, 2);
fs_ma = tbl_freezEpisodes.fs_ma(1);
times_wSpeed = [1:length(wSpeed)]/fs_ma + tdur_trial(1);

% gauss filted
psds_plot = imgaussfilt(psds_AOI,'FilterSize',[3,11]);
clear psds_AOI

% plot spectrogram
T_chnsarea = tbl_freezEpisodes.T_chnsarea{1};
[~, ~, nchns] = size(psds_plot);
for chi = 1: nchns
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
    
    figure();
    
    % plot spectrogram
    subplot(4,1,[2,3,4])
    ax1 = gca;
    imagesc(ax1,times_plot, freqs_plot, psds_plot(:, :, chi));hold on
    
    if ~exist('clim', 'var')|| isempty(clim)
        set(ax1,'YDir','normal')
    else
        set(ax1,'YDir','normal', 'CLim', clim)
    end
    
    colormap(jet)
    colorbar
    
    xlabel('time/s')
    ylabel('Frequency(Hz)')
    
    % change time 0 name
    xtklabels = xticklabels;
    xtklabels{find(cellfun(@(x) strcmp(x,'0'), xtklabels))} = time0name;
    xticklabels(xtklabels);
    plot([0 0], ylim, 'r--', 'LineWidth',1.5)
    pos_ax = get(ax1, 'Position');

    
    % plot wSpeed
    subplot(4,1,1)
    ax2 = gca;
    pos = get(ax2, 'Position');
    pos(3) = pos_ax(3);
    hlegs = [];
    h = plot(times_wSpeed, wSpeed, 'DisplayName', 'wristSpeed');
    hlegs = [hlegs h]; clear h
    hold on
    h = plot(xlim, [30 30], 'k--', 'DisplayName', 'MoveThres');
    hlegs = [hlegs h];
    clear h
    
    set(ax2, 'Position', pos);
    xlim(get(ax1, 'XLim'));
    ylimit = ylim;
    if ylimit(2) < 32
        ylim([0 32])
    end    
    plot([0 0], ylim, 'r--', 'LineWidth',1.5)
    set(gca, 'XTick', get(ax1, 'XTick'), 'XTickLabel', get(ax1, 'XTickLabel'))
    ylabel('speed')
    legend(hlegs, 'Position', [0.82 0.85 0.14 0.07])
    clear pos_ax pos hlegs
    
    
    title(ax2, ['Freeze2Reach-MA-spectrogram-' brainarea ' [' num2str(f_AOI(1)) ' ' num2str(f_AOI(2)) ']Hz'])
    savefile = fullfile(savefolder, ['Freeze2Reach-MA-spectrogram-' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz-' brainarea]);
    print(gcf, savefile, '-painters', '-depsc')
    close all
    clear ax1 ax2
end



end


function [psd_allchns, freqs, times] = calc_psd_allchns(lfp, fs, varargin)
% calculate psd of lfp 
%
%   Input:
%       lfp: nchns * ntemp
%       fs: sample rate
%
%       Name-Value:
%               twin: used in spectrogram, time window for segment (default 0.2 s)
%
%               toverlap: used in spectrogram, time window for overlap (default 0.18 s)
%
%   Return:
%       psd_allchns: nf * nt * nchns
%       freqs: nf * 1
%       times: 1 * nt


% parse 
p = inputParser;
addParameter(p, 'twin',0.5, @(x)isscalar(x)&&isnumeric(x));
addParameter(p, 'toverlap',0.4, @(x)isscalar(x)&&isnumeric(x));
parse(p,varargin{:});
twin = p.Results.twin;
toverlap = p.Results.toverlap;


% calculate psd for each chn 
nwin = round(twin * fs);
noverlap = round(toverlap * fs);
[nchns, ~] = size(lfp);
psd_allchns = [];
for chi = 1 : nchns
    x = lfp(chi, :);
    [~, freqs, times, psd] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
    
    % convert into dB
    psd = 10 * log10(psd);
    
    % append
    psd_allchns = cat(3, psd_allchns, psd); % psd_allchns: nf * nt * nchns
    
    clear x psd
end
end