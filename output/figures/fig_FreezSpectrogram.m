function fig_FreezSpectrogram()
codefilepath = mfilename('fullpath');


% find the codefolder
tmp = regexp(codefilepath, '.*\code', 'match');
if length(tmp) ~= 1
    disp('can not find code path correctly.')
    return;
end
codefolder = tmp{1};
clear tmp

% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));

%% Input & save
[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
inputfolder = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm3_fs500Hz_freezeSKTData_EpisodeExtract');



savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;


%% plot figure parameters
w_colormap = 350; % width  for the colormap
h_colormap = 120; % height for the colormap

w_deltax1_colormap = 5; % x distance between two color map within the same NHP
w_deltax2_colormap = 20; % x distance between two color map of different NHPs

w_textMovePhase = 70; % width showing the moveing phase, i.e. preMove
w_textpair = 80; % width showing the pair name, i.e. M1-STN
w_textColorbar = 80; % width showing the colarbar 

h_deltay_colormap_J = 80; % y distance between two color map of animal J
h_deltay_colormap_K = 10; % y distance between two color map of animal K

h_textAnimal = 30; % height showing the animal name, i.e. animal J/K
h_textCond = 30; % height showing the condition, i.e. Mild-Normal
h_textFreNum = 10; % height showing the frequency number, i.e. 10 12
h_textFreLabel = 40; % height showing the frequency label, i.e. Frequences/Hz


fontname = 'Times New Roman';


%% Code start here

files = dir(fullfile(inputfolder, ['*moderate*.mat']));

% extract freeze tbl_freezEpisodes
[tbl_freezEpisodes]= tbl_freezEpisodes_extract(files, 't_bef', 0.8, 't_aft', 1);

[~, combFreeTypes] = optFreezeTypes_extract();
for cfi = 1: length(combFreeTypes)
    freezeType = combFreeTypes{cfi};
    [psds, freqs_psd, times_psd, wSpeeds, times_wSpeed]= extract_spectrogramMA_acrossTrials(tbl_freezEpisodes, freezeType);
    for chi = 1 : size(psds, 3)
        psds_1chn = squeeze(psds(:, :, chi));
        plot_spectrogramMA_1chn(psds_1chn, freqs_psd, times_psd, wSpeeds, times_wSpeed);
        saveas(gcf, '1.png')
        clear psds_1chn
    end
    
    clear freezeType
end



function [psds, freqs_psd, times_psd, wSpeeds, times_wSpeed]= extract_spectrogramMA_acrossTrials(tbl_freezEpisodes, freezeType)
%
%
%  Return:
%       psds: nfs * nts_psd * nchns 
%       freqs_psd: frequencies nfs * 1
%       times_psd: 1 * nts_psd
%       wSpeeds: speed of all trials nts_speed * ntrials
%       times_wSpeed : time point for speed 1 * nts_speed

tdur = 5;
f_AOI = [8 40];

tbl_subfreezEpi = tbl_freezEpisodes(tbl_freezEpisodes.freezType == freezeType, :);

% calc each freeze phase
psds_alltrials = []; % nf * nt * nchns * ntrials
wSpeeds = []; % wSpeeds: ntemp * ntrials
for tbi = 1 : height(tbl_subfreezEpi)
    
    % extract lfp_freeze duration
    lfp = tbl_subfreezEpi.lfp{tbi}; % lfp: nchns * ntemp
    fs_lfp = tbl_subfreezEpi.fs_lfp(tbi);

    idx_lfpStrendFreez = tbl_subfreezEpi.idx_lfpStrendFreez{tbi};
    
    t_freeze = (idx_lfpStrendFreez(2) - idx_lfpStrendFreez(1))/fs_lfp;
    if t_freeze < tdur
        disp(['less than ' num2str(tdur)])
        continue;
    end
    lfp = lfp(:,1: round(idx_lfpStrendFreez(1)+fs_lfp*tdur));
    
    % extract wSpeed_freeze duration
    fs_ma = tbl_subfreezEpi.fs_ma(tbi);
    idx_wSpeedStrendFreez = tbl_subfreezEpi.idx_wSpeedStrendFreez{tbi};
    wSpeed = tbl_subfreezEpi.wSpeed{tbi};
    wSpeed = wSpeed(1: round(idx_wSpeedStrendFreez(1)+fs_ma*tdur), 1); % wSpeed: ntemp * 1
    wSpeeds = cat(2, wSpeeds, wSpeed);
    clear fs_ma  idx_wSpeedStrendFreez wSpeed
    
    % calc psd
    [psd, freqs, times] = calc_psd_allchns(lfp, fs_lfp);% psds: nf * nt * nchns
    
    if tbi == 1
        % extract based on f_AOI and align times with t_lfpfreezstr as time 0s
        idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        freqs_psd =  freqs(idx_f);
        idx_lfpfreezstr = idx_lfpStrendFreez(1);
        times_psd = times - idx_lfpfreezstr/fs_lfp;
        clear idx_lfpfreezstr   
    end
    psd = psd(idx_f, :, :);
    
    % append
    psds_alltrials = cat(4,psds_alltrials,psd);
   
    clear psds lfp idx_lfpStrendFreez psd freqs times fs_lfp
end
psds = mean(psds_alltrials, 4);


fs_ma = tbl_freezEpisodes.fs_ma(1);
idx_wSpeedfreezstr = tbl_subfreezEpi.idx_wSpeedStrendFreez{tbi}(1);
times_wSpeed = [1:size(wSpeeds, 1)]/fs_ma - idx_wSpeedfreezstr/fs_ma;



function plot_spectrogramMA_1chn(psds_1chn, freqs_psd, times_psd, wSpeeds, times_wSpeed, varargin)
%
%
%  Inputs:
%       psds: nfs * nts_psd 
%       freqs_psd: frequencies nfs * 1
%       times_psd: 1 * nts_psd
%       wSpeeds: speed of all trials nts_speed * ntrials
%       times_wSpeed : time point for speed 1 * nts_speed
%
%       Name-Value: 
%           'clim' - clim (default [])
%           'show_xticklabels_spect' - show (true, default) or not show (false) xticklabels for spectrogram
%           'show_xlabel_spect' - show (true, default) or not show (false) xlabel for spectrogram
%           'show_ylabel_spect' - show (true, default) or not show (false) ylabel for spectrogram
%           'show_colorbar_spect' - show (true, default) or not show (false) colorbar for spectrogram
%           'show_titlename_speed' - show (true, default) or not show (false) titlename for speed
%           'show_ylabel_speed' - show (true, default) or not show (false) ylabel for speed
%           'show_MoveThres_speed' - show (true, default) or not show (false) MoveThres line for speed
%           'show_legends_speed' - show (true, default) or not show (false) legend for speed
%           'time0name' - text name for time 0, default '0'
%           'titlename_speed' - title name for speed, default 'Speed'



% parse params
p = inputParser;
addParameter(p, 'clim', [], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'show_xticklabels_spect', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xlabel_spect', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel_spect', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar_spect', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_titlename_speed', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel_speed', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_MoveThres_speed', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_legends_speed', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'time0name', '0', @isstr);
addParameter(p, 'titlename_speed', 'Speed', @isstr);


parse(p,varargin{:});
show_xticklabels_spect = p.Results.show_xticklabels_spect;
show_xlabel_spect = p.Results.show_xlabel_spect;
show_ylabel_spect = p.Results.show_ylabel_spect;
show_colorbar_spect = p.Results.show_colorbar_spect;
show_titlename_speed = p.Results.show_titlename_speed;
show_ylabel_speed = p.Results.show_ylabel_speed;
show_MoveThres_speed = p.Results.show_MoveThres_speed;
show_legends_speed = p.Results.show_legends_speed;
time0name = p.Results.time0name;
titlename_speed = p.Results.titlename_speed;




%%% plot 
fig = figure();

%%% subplot spectrogram

%gauss filted
psds_1chn = imgaussfilt(psds_1chn,'FilterSize',[3,11]);

ax_spec = subplot(4,1,[2,3,4]);
imagesc(ax_spec,times_psd, freqs_psd, psds_1chn);hold on

if ~exist('clim', 'var')|| isempty(clim)
    set(ax_spec,'YDir','normal')
else
    set(ax_spec,'YDir','normal', 'CLim', clim)
end
colormap(jet)
colorbar


% plot time 0 line
plot(ax_spec, [0 0], ylim, 'k--')

pos_spec = get(ax_spec, 'Position');

% show inf
if show_xticklabels_spect
    % change time 0 name
    xtklabels = xticklabels(ax_spec);
    xtklabels{find(cellfun(@(x) strcmp(x,'0'), xtklabels))} = time0name;
    xticklabels(ax_spec, xtklabels);
else
    xticks([]);
end


if show_xlabel_spect
    xlabel(ax_spec, 'Time (s)', 'fontsize', 12, 'FontName', 'Times New Roma', 'FontWeight', 'bold')
end
if show_ylabel_spect
    ylabel(ax_spec, 'Frequency (Hz)', 'fontsize', 12, 'FontName', 'Times New Roma', 'FontWeight', 'bold')
end

if show_colorbar_spect
    colorbar;
    set(ax_spec, 'Position', pos_spec);
    clear pos
end

%%% subplot wSpeed
ax_speed = subplot(4,1,1);
hlegs = [];
for tri = 1 : size(wSpeeds, 2)
    wSpeed = wSpeeds(:, tri);
    if tri == 1
        h = plot(ax_speed, times_wSpeed, wSpeed, 'DisplayName', 'wristSpeed');
        hold on
        hlegs = [hlegs h];
        clear h
    else
        plot(ax_speed, times_wSpeed, wSpeed);
        hold on
    end
    clear wSpeed
end

% time 0 line
plot([0 0], ylim, 'k--')

% adjust pos, xlim, ylim
pos_speed = get(ax_speed, 'Position');
pos_speed(3) = pos_spec(3);
set(ax_speed, 'Position', pos_speed)
xlim(get(ax_spec, 'XLim'));
ylimit = ylim;
if ylimit(2) < 32
    ylim([0 32])
end

set(ax_speed, 'XTick', [])


% show inf speed
if show_titlename_speed
    title(ax_speed, titlename_speed, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roma')
end

if show_ylabel_speed
    ylabel(ax_speed, 'Speed', 'fontsize', 12, 'FontName', 'Times New Roma', 'FontWeight', 'bold')
end

if show_MoveThres_speed
    h = plot(ax_speed, xlim, [30 30], 'r--', 'DisplayName', 'MoveThres');
    hlegs = [hlegs h];
    clear h
end

if show_legends_speed
    legend(ax_speed, hlegs, 'Position', [0.82 0.85 0.14 0.07])
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
        idx_wSpeedfreezeInExtract = [n_bef+1 length(lfpfreez)-n_aft];
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