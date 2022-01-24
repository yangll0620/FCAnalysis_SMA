function m4_restData_PSD()
%%   PSD estimates for mild and normal individually
%
%       psd for each brain area, as well as each DBS contact
%
%

%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code') - 1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder, 'util')));

% the corresponding pipeline folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%%  input setup

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

% pwelch psd estimate variable
twin_pwelch = 2;

% variables for plotting
plotF_AOI = [5 50];

% input folder
inputfolder = fullfile(codecorresParentfolder, 'm3_restData_rmChns_avgArea');

%% save setup
savefolder = codecorresfolder;
savefilename_prefix = 'psd_';

%% Code Start Here

file_psdall = fullfile(savefolder, [savefilename_prefix '_allsegs_normalmildmoderate.mat']);

%%%  calculate/load dbs psd for normal, mild and moderate %%%
if ~exist(file_psdall, 'file')% not exist
    
    [pxxs_allfiles_normal, F_pxx_normal] = pxx_eacharea_allfiles(dir(fullfile(inputfolder, '*_normal_*.mat')), twin_pwelch);
    [pxxs_allfiles_mild, F_pxx_mild] = pxx_eacharea_allfiles(dir(fullfile(inputfolder, '*_mild_*.mat')), twin_pwelch);
    [pxxs_allfiles_moderate, F_pxx_moderate] = pxx_eacharea_allfiles(dir(fullfile(inputfolder, '*_moderate_*.mat')), twin_pwelch);
    
    if ~isequal(F_pxx_normal, F_pxx_mild, F_pxx_moderate)
        disp('F_pxx_normal, F_pxx_mild and F_pxx_moderate not equal');
        
        return;
    else
        F_pxx = F_pxx_normal;
        
        save(file_psdall, 'pxxs_allfiles_normal', 'pxxs_allfiles_mild', 'pxxs_allfiles_moderate', 'F_pxx');
        clear F_pxx_normal F_pxx_mild F_pxx_moderate
    end
    
else
    load(file_psdall, 'pxxs_allfiles_normal', 'pxxs_allfiles_mild', 'pxxs_allfiles_moderate', 'F_pxx');
end

%%%  plot  %%%
brainareas = fieldnames(pxxs_allfiles_mild);

for i = 1:length(brainareas)
    brainarea = brainareas{i};
    
    % load normal and mild data
    eval(['psd_normal = pxxs_allfiles_normal.' brainarea ';'])
    eval(['psd_mild = pxxs_allfiles_mild.' brainarea ';'])
    eval(['psd_moderate = pxxs_allfiles_moderate.' brainarea ';'])
    
    %  psd_normal, psd_mild: nfs * nsegs
    plotPSD_comp_1chn(psd_normal, psd_mild, psd_moderate, F_pxx, plotF_AOI, savefolder, brainarea, animal)
end

%%% combine %%%
close all

end

function [pxxs_allfiles, F_pxx] = pxx_eacharea_allfiles(files, twin_pwelch)
%% extract psd from all the files, each psd for each dbs contact and one psd for one area (except dbs)
%
% Arg:
%       files = dir(fullfile('.', '*_mild_*.mat'));
%       twin_pwelch : twin for pwelch
%
% Outputs:
%
%       pxxs_allfiles: the PSD estimate of all the segments from the files
%           e.g. pxxs =
%                   struct with fields:
%                      M1: [nfs * nsegs double]
%                     stn0_1: [nfs * nsegs double]
%                      gp_0-1: [nfs * nsegs double]
%
%       F_pxx: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)

nfiles = length(files);

for filei = 1:nfiles
    
    file = fullfile(files(filei).folder, files(filei).name);
    
    [pxxs_1file, F_pxx_1file] = pxx_eacharea_onefile(file, twin_pwelch);
    
    brainareas_saved = fieldnames(pxxs_1file);
    
    % combined pxxs from all the files
    if (~exist('pxxs_allfiles', 'var'))
        pxxs_allfiles = pxxs_1file;
    else
        
        for i = 1:length(brainareas_saved)
            brainarea = brainareas_saved{i};
            %  pxxs_allfiles.M1, pxxs.M1: nfs * nsegs
            eval(['pxxs_allfiles.' brainarea '= cat(2, pxxs_allfiles.' brainarea ', pxxs_1file.' brainarea ');'])
            
            clear onearea
        end
        
    end
    
    % extract the F_pxx
    if ~exist('F_pxx', 'var')
        F_pxx = F_pxx_1file;
    else
        
        if F_pxx ~= F_pxx_1file
            disp(['F_pxx ~=f for ' brainarea ', segi = ' num2str(segi) '']);
            
            F_pxx = [];
            pxxs_allfiles = [];
            
            return;
        end
        
    end
    
    clear file pxxs_1file F_pxx_1file brainareas
end

end

function [pxxs, F_pxx] = pxx_eacharea_onefile(file, twin_pwelch)
%% extract psd of all the segments from the files, each psd for each dbs contact and one psd for one area (except dbs)
%
% Outputs:
%
%       pxxs: the PSD estimate of all the segments from the files
%           e.g. pxxs =
%                   struct with fields:
%                      M1: [nfs * nsegs double]
%                     stn0_1: [nfs * nsegs double]
%                      gp0_1: [nfs * nsegs double]
%
%       F_pxx: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)

% load data
load(file, 'fs', 'data_segments', 'T_chnsarea');

if (isempty(data_segments))% data_segments is empty
    F_pxx = [];
    pxxs = [];
end



% psd pwelch paramers
nwins = round(twin_pwelch * fs);
noverlap = round(nwins * 0.9);

pxxs = struct();

for chi = 1:height(T_chnsarea)
    
    brainarea = T_chnsarea.brainarea{chi};
    brainarea_saved = strrep(brainarea,'-','_');% change stn0-1 to stn0_1
    for segi = 1:length(data_segments)
        
        % extract the lfp data of segi , lfp_oneseg: 1* ntemp
        lfp_oneseg = data_segments(segi).lfp(chi, :);
        
        
        %%% calcualte PSD throuch zscore and pwelch %%%
        
        % zscore of the lfp_oneseg along each column
        lfp_zscore = zscore(lfp_oneseg);
        
        % PSD is computed independently for each column, Pxx: nfs * nchns, f: nfs * 1
        [Pxx, f] = pwelch(lfp_zscore, nwins, noverlap, nwins, fs);
        
        if ~exist('F_pxx', 'var')
            F_pxx = f;
        else
            
            if F_pxx ~= f
                disp(['F_pxx ~=f for ' brainarea ', segi = ' num2str(segi) '']);
                
                F_pxx = [];
                pxxs = [];
                
                return;
            end
            
        end
        
        
        if ~isfield(pxxs, brainarea_saved)
            % first time calculate pxx for brainarea
            eval(['pxxs.' brainarea_saved ' = Pxx;'])
        else
            % Pxx: nfs * 1; pxxs.M1: nfs * nsegs
            eval(['pxxs.' brainarea_saved ' = cat(2, pxxs.' brainarea_saved ', Pxx);'])
        end
        
        
        clear lfp_oneseg lfp_zscore Pxx f
    end
    
    clear brainarea mask_area segi brainarea_saved
end

end

function plotPSD_comp_1chn(psd_normal, psd_mild, psd_moderate, F_all, plotF_AOI, savefolder, brainarea, animal)
%%  plot the psd comparison of normal and mild
%
%   Inputs:
%       psd_normal, psd_mild: psd of all segments in normal or mild, nfs * nsegs
%
%       F_all: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
%
%       plotF_AOI: the plot frequence range (in hertz)  (1 * 2)
%
%       savefile_prefix: the savefile prefix including folder

% find the idx for F_AOI
idx_AOI = find(F_all >= plotF_AOI(1) & F_all <= plotF_AOI(2));

% colors setup
color_normal_range = [224, 255, 255] / 255;
color_normal_mean = [0, 0, 255] / 255;
color_mild_range = [255, 228, 225] / 255;
color_mild_mean = [255, 00, 0] / 255;
color_moderate_range = [238, 238, 238] / 255;
color_moderate_mean = [0, 0, 0] / 255;

% plot setup
linewidth = 1.5;

% F_AOI
F_AOI = F_all(idx_AOI);
% reshape F_AOI
[m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m * n); clear m n

% extract the  psd of AOI
psd_normal_FAOI = psd_normal(idx_AOI, :);
psd_mild_FAOI = psd_mild(idx_AOI, :);
psd_moderate_FAOI = psd_moderate(idx_AOI, :);

psd_normal_high = max(psd_normal_FAOI, [], 2);
psd_normal_low = min(psd_normal_FAOI, [], 2);
psd_normal_mean = mean(psd_normal_FAOI, 2);

psd_mild_high = max(psd_mild_FAOI, [], 2);
psd_mild_low = min(psd_mild_FAOI, [], 2);
psd_mild_mean = mean(psd_mild_FAOI, 2);


psd_moderate_high = max(psd_moderate_FAOI, [], 2);
psd_moderate_low = min(psd_moderate_FAOI, [], 2);
psd_moderate_mean = mean(psd_moderate_FAOI, 2);

% reshape, psd_*_high/low/mean into 1 * nfs

[m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m * n); clear m n
[m, n] = size(psd_normal_low); psd_normal_low = reshape(psd_normal_low, 1, m * n); clear m n
[m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m * n); clear m n

[m, n] = size(psd_mild_high); psd_mild_high = reshape(psd_mild_high, 1, m * n); clear m n
[m, n] = size(psd_mild_low); psd_mild_low = reshape(psd_mild_low, 1, m * n); clear m n
[m, n] = size(psd_mild_mean); psd_mild_mean = reshape(psd_mild_mean, 1, m * n); clear m n


[m, n] = size(psd_moderate_high); psd_moderate_high = reshape(psd_moderate_high, 1, m * n); clear m n
[m, n] = size(psd_moderate_low); psd_moderate_low = reshape(psd_moderate_low, 1, m * n); clear m n
[m, n] = size(psd_moderate_mean); psd_moderate_mean = reshape(psd_moderate_mean, 1, m * n); clear m n

% plot range
figure
set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
fill([F_AOI flip(F_AOI)], [psd_normal_high flip(psd_normal_low)], color_normal_range, 'LineStyle', 'none')
hold all
fill([F_AOI flip(F_AOI)], [psd_mild_high flip(psd_mild_low)], color_mild_range, 'LineStyle', 'none')
fill([F_AOI flip(F_AOI)], [psd_moderate_high flip(psd_moderate_low)], color_moderate_range, 'LineStyle', 'none')

% plot mean
h1 = plot(F_AOI, psd_normal_mean, 'Color', color_normal_mean, 'LineWidth', linewidth);
h2 = plot(F_AOI, psd_mild_mean, 'Color', color_mild_mean, 'LineWidth', linewidth);
h3 = plot(F_AOI, psd_moderate_mean, 'Color', color_moderate_mean, 'LineWidth', linewidth);
%
%     % find the frequency with maximum density
%     [maxPSD, idx_max] = max(psd_mild_mean);
%     F_maxPSD = round(F_AOI(idx_max));
%     plot([F_maxPSD F_maxPSD], [0 maxPSD + maxPSD * 0.2], 'k--')

xlim([min(F_AOI) max(F_AOI)])
ylim([0 0.2])
set(gca,'XTick',[10 15 20 25 30 35 40 45 50],'YTick',[0 0.1 0.2]);
xlabel('Frequency (Hz)', 'FontWeight','bold')
ylabel('Power', 'FontWeight','bold')
set(gca, 'Box', 'off')

% legend
legend([h1, h2, h3], {'normal', 'mild', 'moderate'})

% title
title([animal ' Rest PSD in ' strrep(brainarea, '_', '-')])

% save figure
savename = fullfile(savefolder, [animal 'Rest_psd_' brainarea]);
saveas(gcf, savename, 'png')

clear psd_allsegs_normal psd_allsegs_mild psd_allsegs_moderate
clear psd_normal_FAOI psd_mild_FAOI psd_moderate_FAOI
clear psd_normal_high psd_normal_low psd_normal_mean psd_mild_high psd_mild_low psd_mild_mean psd_moderate_high psd_moderate_low psd_moderate_mean
clear h1 h2 h3 maxPSD F_maxPSD idx_max
clear savename
end

