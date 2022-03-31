function [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= freez_lfpsegs(files, varargin)
% extract freeze lfpsegs
%   Inputs
%       Name-Value: 
%           't_ThreFreeze': time threshold for extracted freezing phase,
%                           only extract the segs longer than the threshold (default = 5)
%
%           'tdur': extracted seg duration, default = [0 t_ThreFreeze], 0
%                   is the freeze time onset

if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpsegs_freeze = [];
    fs = [];
    T_chnsarea = [];
    combFreeTypes = [];
    
    return;
end


p = inputParser;
addParameter(p, 'tdur', [], @(x)isempty(x)||(isnumeric(x)&&isvector(x)&&length(x)==2));
addParameter(p, 't_ThreFreeze',5, @(x)isscalar(x)&&isnumeric(x));
parse(p,varargin{:});
tdur = p.Results.tdur;
t_ThreFreeze = p.Results.t_ThreFreeze;
if isempty(tdur)
    tdur = [0 t_ThreFreeze];
end


optFreezeTypes = optFreezeTypes_extract();
combFreeTypes = {'InitFreeze', 'ReachFreeze', 'ManipuFreeze'}; % combined {'freeze during React-Reach'}  and  {'freeze during Reach'} 

lfpsegs_freeze = struct();
for frTi = 1 : length(combFreeTypes)
    lfpsegs_freeze.(combFreeTypes{frTi}) = [];
end

for fi = 1: length(files)
    loadfilename = files(fi).name;
    load(fullfile(files(fi).folder, loadfilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'freezStruct', 'selectedTrials');
    
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
        
        t_str = freezEpisodes{frzi}.freezeTPhaseS(1);
        t_end = freezEpisodes{frzi}.freezeTPhaseS(2);
        if t_end - t_str < t_ThreFreeze
            clear tri t_str t_end
            continue;
        end
        idx_dur = round((tdur + t_str) * fs_lfp);
        if idx_dur(1) == 0
            idx_dur(1) = 1;
        end
        seglfp  = lfpdata{tri}(:, idx_dur(1): idx_dur(2));

        
        %%% append to lfpsegs_freeze
        freezeType = freezEpisodes{frzi}.freezeType;
        frTi = find(strcmp(freezeType, optFreezeTypes));
        if frTi ==3 || frTi == 4
            frTi = frTi -1;
        end
        lfpsegs_freeze.(combFreeTypes{frTi}) = cat(3, lfpsegs_freeze.(combFreeTypes{frTi}), seglfp);
        
        clear tri t_str t_end freezeType frTi seglfp
        clear idx_dur
    end   
    clear freezEpisodes
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;
end


function plot_Spectrogram_acrossSegs(lfpsegs, T_chnsarea, tdur, fs, varargin)
% plot lfpdata of all the channels: nchns * ntemp * nsegs




twin = 0.2;
toverlap = 0.18;

nwin = round(twin * fs);
noverlap = round(toverlap * fs);


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_AOI', [-0.5 0.5], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 'clim', [], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 'title_prefix', '', @isstr);
addParameter(p, 'savefilename_prefix', '', @isstr);
addParameter(p, 'savefolder', '', @isstr);
addParameter(p, 'img_format', 'tif', @isstr);

parse(p,varargin{:});
f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
clim = p.Results.clim;
img_format = p.Results.img_format;
title_prefix = p.Results.title_prefix;
savefilename_prefix = p.Results.savefilename_prefix;
savefolder = p.Results.savefolder;
if isempty(title_prefix)
    title_prefix = 'spectrogram';
end
if isempty(savefilename_prefix)
    savefilename_prefix = 'spectrogram';
end
if isempty(savefolder)
    savefolder = pwd;
end


% calculate psd for each chn across trials
[nchns, ~, nsegs] = size(lfpsegs);
psd_allchns = [];
for chi = 1 : nchns
    psds = []; %  psds: nf * nt * ntrials
    for segi = 1: nsegs
        x = lfpsegs(chi, :, segi);
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
    
    times = times + tdur(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end

% plot spectrogram
[nf, nt, nchns] = size(psd_allchns);
for chi = 1: nchns
    brainarea = T_chnsarea.brainarea{chi};
    if strcmp(brainarea, 'M1')
        clim = [-35 -10];
    end
    if contains(brainarea, 'stn')
        clim = [-30 -10];
    end
    if contains(brainarea, 'gp')
        clim = [-35 -15];
    end
    
    figure();
    ax = gca;
    imagesc(ax,times_plot, freqs_plot, psd_allchns(:, :, chi));
    
    if ~exist('clim', 'var')|| isempty(clim)
        set(ax,'YDir','normal')
    else
        set(ax,'YDir','normal', 'CLim', clim)
    end
    
    colormap(jet)
    colorbar
    
    xlabel('time/s')
    ylabel('Frequency(Hz)')
    
    title(ax, [title_prefix '-' brainarea, ', nsegs=' num2str(nsegs)])
    savefile = fullfile(savefolder, [savefilename_prefix '_' brainarea]);
    saveas(gcf, savefile, img_format);
    close all
end

end