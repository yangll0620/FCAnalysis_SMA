function plot_spectrogram_acrossTrials(tbl_freezEpisodes, freezeType, f_AOI, savefolder)
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

img_format = 'tif';
tbl_subfreezEpi = tbl_freezEpisodes(tbl_freezEpisodes.freezType == freezeType, :);

switch freezeType
    case 'initFreeze'
        time0name = 'cueonset';
        timeEndname = 'freezeEnd';
    case 'reachFreeze'
        time0name = 'freezeStart';
        timeEndname = 'freezeEnd';
    case 'maniFreeze'
        time0name = 'touch';
        timeEndname = 'mouth';
end

tdur = 5;

% calc each freeze phase
psds_alltrials = [];
for tbi = 1 : height(tbl_subfreezEpi)
    lfp = tbl_subfreezEpi.lfp{tbi}; % lfp: nchns * ntemp
    fs = tbl_subfreezEpi.fs_lfp(tbi);
    T_chnsarea = tbl_subfreezEpi.T_chnsarea{tbi};

    idx_strendFreez = tbl_subfreezEpi.idx_strendFreez{tbi};
    
    t_freeze = (idx_strendFreez(2) - idx_strendFreez(1))/fs;
    if t_freeze < tdur
        disp(['less than ' num2str(tdur)])
        continue;
    end
    lfp = lfp(:,1: round(idx_strendFreez(1)+fs*tdur));
    
    % calc psd
    [psds, freqs, times] = calc_psd_allchns(lfp, fs);% psds: nf * nt * nchns
    
    % append
    psds_alltrials = cat(4,psds_alltrials,psds);
   
    clear psds lfp
end
psds = mean(psds_alltrials, 4);
    
% extract based on f_AOI and align times with t_freezstr as time 0s
idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
freqs_plot =  freqs(idx_f);
psd_AOI = psds(idx_f, :, :);
idx_freezstr = idx_strendFreez(1);
times_plot = times - idx_freezstr/fs;
clear idx_f idx_freezstr
    
% gauss filted
psd_plot = imgaussfilt(psd_AOI,'FilterSize',[3,11]);
clear psd_AOI
    
% plot spectrogram
title_prefix = [freezeType '-acrossTrials'];
[~, ~, nchns] = size(psd_plot);
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
    imagesc(ax,times_plot, freqs_plot, psd_plot(:, :, chi));hold on
    
    if ~exist('clim', 'var')|| isempty(clim)
        set(ax,'YDir','normal')
    else
        set(ax,'YDir','normal', 'CLim', clim)
    end
    
    colormap(jet)
    colorbar
    
    xlabel('time/s')
    ylabel('Frequency(Hz)')
    
    % change time 0 name
    xtklabels = xticklabels;
    xtklabels{find(cellfun(@(x) strcmp(x,'0'), xtklabels))} = time0name;
    xticklabels(xtklabels);
    plot([0 0], ylim, 'k--')
    
    title(ax, [title_prefix '-' brainarea])
    savefile = fullfile(savefolder, [title_prefix '_' brainarea]);
    saveas(gcf, savefile, img_format);
    close all
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
%               toverlap: used in spectrogram, time window for overlap (default 0.18 s)
%
%   Return:
%       psd_allchns: nf * nt * nchns
%       freqs: nf * 1
%       times: 1 * nt


% parse 
p = inputParser;
addParameter(p, 'twin',0.2, @(x)isscalar(x)&&isnumeric(x));
addParameter(p, 'toverlap',0.18, @(x)isscalar(x)&&isnumeric(x));
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