function m2_chnOfI_restData_plotSpectrogram15s()
%%   plot 15s spectrogram for moderate Kitty Rest data for chns of Interest


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
inputfolder = fullfile(codecorresParentfolder, 'm1_restData_rmChns_avgArea');

pdcond = 'moderate';

% paras for calc and plot spectrogram
twin = 1;
toverlap = 0.8;
f_AOI = [8 40];

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

saveimg_format = 'jpg';

%% Code Start Here
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
for fi = 1 : length(files) 
    filename = files(fi).name;
    load(fullfile(inputfolder, filename), 'data_segments', 'fs', 'T_chnsarea');
    
    % extract datebkstr
    datebkstr = regexp(filename, '\d{8}_tdt\d+', 'match');
    datebkstr = strrep(datebkstr{1}, '_', '-');
    
    % extract mask_UsedChns
    mask_usedChns = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    
    for segi = 1 : length(data_segments)
        lfp = data_segments(segi).lfp;
        lfp = lfp(mask_usedChns, :);
        
        [psds_allchns, freqs, times] = calc_psd_allchns(lfp, fs, 'twin', twin, 'toverlap', toverlap);
        
        % select freqs and corresponding psd
        idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        freqs_AOI =  freqs(idx_f);
        psds_allchns_AOI = psds_allchns(idx_f, :, :);
        times_AOI = times;
        
        % plot and save
        for chi = 1 : size(psds_allchns_AOI, 3)
            psds_plot = psds_allchns_AOI(:, :, chi);
            
            % gaussian filter
            psds_plot = imgaussfilt(psds_plot,1, 'FilterSize',5);
            
            brainarea = T_chnsarea.brainarea{chi};
            if contains(brainarea, 'gp')
                clim = [-160 -100];
            end
            if strcmp(brainarea, 'M1')
                clim = [-140 -85];
            end
            if contains(brainarea, 'stn')
                clim = [-160 -95];
            end

            
            
            savefilename = [animal '-' pdcond '-' brainarea '-' datebkstr  '-seg' num2str(segi) '-'];
            titlename = [animal '-' pdcond '-' brainarea '-' datebkstr  '-seg' num2str(segi)];
            
            plotsave_spectrogram(psds_plot, freqs_AOI, times_AOI, savefolder, savefilename, titlename, 'clim', clim,'savedimg_format', 'tif',  'time0name', '');
            
            clear psds_plot 
        end
        
        
        clear lfp psd_allchns freqs times
        clear idx_f  freqs_AOI psds_allchns_AOI times_AOI
    end
    
    % clear
    clear mask_usedChns T_chnsarea 
    clear('data_segments', 'fs', 'T_chnsarea');
    clear filename datebkstr
end
end


function plotsave_spectrogram(psds, freqs, times, savefolder, savefilename, titlename, varargin)
% plot one channel psds with freqs and times  
%
%   Usage:
%       plotsave_spectrogram(psds, freqs, times, savefolder, savefilename, titlename, 'savedimg_format', 'tif', 'clim', [-35 -10], 'time0name', '');
%
%   Input:
%       psds: nfreqs * ntimes
%       freqs: corresponding frequences nfreqs * 1
%       times: corresponding times 1 * ntimes
%       savefolder: the saved into folder
%
%       Name-Value:
%               savedimg_format: the saved image format, default 'tif'
%
%               clim: used in spectrogram, default []
%
%   Save:
%       save with savefilename in savefolder


% parse 
p = inputParser;
addParameter(p, 'savedimg_format','tif', @(x)ischar(x));
addParameter(p, 'clim',[], @(x)isempty(x)||(isnumeric(x)&&isvector(x)&&length(x)==2));
addParameter(p, 'time0name','', @(x)ischar(x));
parse(p,varargin{:});
savedimg_format = p.Results.savedimg_format;
clim = p.Results.clim;
time0name = p.Results.time0name;

if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end


% plot
figure();
ax = gca;
imagesc(ax,times, freqs, psds);
hold on

if ~exist('clim', 'var')|| isempty(clim)
    set(ax,'YDir','normal')
else
    set(ax,'YDir','normal', 'CLim', clim)
end

colormap(jet)
colorbar

xlabel('time/s')
ylabel('Frequency(Hz)')

if ~isempty(time0name)
    % change time 0 name
    xtklabels = xticklabels;
    xtklabels{find(cellfun(@(x) strcmp(x,'0'), xtklabels))} = time0name;
    xticklabels(xtklabels);
    plot([0 0], ylim, 'k--')
end



title(ax, titlename)
savefile = fullfile(savefolder, savefilename);
saveas(gcf, savefile, savedimg_format);
close all

end



function [psd_allchns, freqs, times] = calc_psd_allchns(lfp, fs, varargin)
% calculate psd of lfp 
%
%   Usage:
%       [psd_allchns, freqs, times] = calc_psd_allchns(lfp, fs, 'twin', 1, 'toverlap', 0.8);
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
