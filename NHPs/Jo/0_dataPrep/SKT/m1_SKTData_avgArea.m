function m1_SKTData_avgArea()
%% average lfp across each area
%
%
%   1. remove chns marked with multiple areas or GP/STN in Gray Matter
%   2. average across each area (e.g M1)
%   3. change STN/GP into stn0-1, stn1-2... gp0-1, gp1-2 et.al



codefilepath = mfilename('fullpath');

% codefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path and NHP path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% datafolder, pipelinefolder 
[correspipelinefolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


% animal
animal = animal_extract(correspipelinefolder);


%%  input setup

% input folder: extracted raw STK data 
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');


%% save setup
savefolder = correspipelinefolder;
savecodefolder = fullfile(correspipelinefolder, 'code');
savefilename_addstr = 'avgArea';

%% start here
if ~exist(savecodefolder, 'dir')
    mkdir(savecodefolder);
end
copyfile2folder(codefilepath, savecodefolder);


files = dir(fullfile(inputfolder, '*.mat'));


nfiles = length(files);
close all;
for fi = 22 : nfiles

    disp(['Avg area lfp data in file ' num2str(fi) '/' num2str(nfiles)]);

    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(fi).name;
    load(fullfile(files(fi).folder, filename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'T_idxevent_lfp', ...
                                               'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    
    if ~exist('lfpdata', 'var')
        continue;
    end
                                           
    % remove chns marked with multiple areas or GP/STN in Gray Matter
    multipleAreas_mask = cellfun(@(x) contains(x, '/'), T_chnsarea.brainarea);
    GPSTNinGM_mask = cellfun(@(x) contains(x, 'GP'), T_chnsarea.brainarea) & ~strcmp(T_chnsarea.electype, 'DBS');
    T_chnsarea = T_chnsarea(~multipleAreas_mask & ~GPSTNinGM_mask, :); T_chnsarea.chni = [1:height(T_chnsarea)]';
    lfpdata = lfpdata(~multipleAreas_mask, :, :);
    
    
    % average across each GM area
    mask_notDBS = ~strcmp(T_chnsarea.electype, 'DBS');
    mask_DBS = strcmp(T_chnsarea.electype, 'DBS');
    lfpdata_GM = lfpdata(mask_notDBS, :, :);
    lfpdata_DBS = lfpdata(mask_DBS, :, :);
    T_notDBSchnsarea = T_chnsarea(mask_notDBS, :);
    T_DBSchnsarea = T_chnsarea(mask_DBS, :);
    [avglfp_notDBS, T_notDBSchnsarea_new] = avglfp_acrossArea(lfpdata_GM, T_notDBSchnsarea);
    
    if isempty(avglfp_notDBS)
        continue;
    end
    
    
    % change STN/GP into stn0-1, stn1-2... gp0-1, gp1-2 et.al
    for i = 1: height(T_DBSchnsarea)
        chn = T_DBSchnsarea(i, :).recordingchn;
        barea = cell2mat(T_DBSchnsarea(i, :).brainarea);
        T_DBSchnsarea(i, :).brainarea = {[lower(barea) num2str(chn-1) '-' num2str(chn)]};
        clear chn
    end
    
    
    % cat GM and DBS
    lfpdata = cat(1, avglfp_notDBS, lfpdata_DBS);
    T_chnsarea = [T_notDBSchnsarea_new; T_DBSchnsarea];
    T_chnsarea.chni = [1: height(T_chnsarea)]';
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx + tmpn - 1) savefilename_addstr ...
        upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];
    
    save(fullfile(savefolder, savefilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'T_idxevent_lfp', ...
                                               'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    
    
    
    clear lfpdata fs_lfp T_chnsarea T_idxevent_lfp fs_ma  T_idxevent_ma smoothWspeed_trial Wpos_smooth_trial Wrist_smooth_trial
    clear mask_GM mask_DBS lfpdata_GM lfpdata_DBS T_GMchnsarea T_DBSchnsarea
    clear avglfp_GM T_GMchnsarea_new
end
end

function [avglfp, T_chnsarea_new] = avglfp_acrossArea(lfpdata, T_GMchnsarea)
% average lfp across each area 
% 
% args:
%       lfpdata : lfp data in area (nchns * ntemp * ntrials )
%       T_chnsarea: chns-area table 
% return:
%       avglfp: averaged lfp data (nareas * ntemp * ntrials
%       T_chnsarea_new: the new chns-area table


uniqAreas = unique(T_GMchnsarea.brainarea);
avglfp = [];
T_chnsarea_new = T_GMchnsarea([], :);
for ai = 1 : length(uniqAreas)
    GMarea = uniqAreas{ai};
    
    
    mask_area = strcmp(T_GMchnsarea.brainarea, GMarea);
    lfp_area = lfpdata(mask_area, :, :);
    avglfp_area = mean(lfp_area, 1);
    
    
    avglfp = cat(1, avglfp, avglfp_area);
    
    T_chnsarea_new = [T_chnsarea_new; {ai, GMarea, nan,  'Gray Matter', []}];
    
    clear GMarea mask_area avglfp_area
end



end

