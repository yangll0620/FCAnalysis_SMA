function m1_SKTData_avgArea()
%% average lfp across each area
%
%   remove unwanted stn3-4, stn4-5, stn5-6 and stn6-7
%
%   averaged across M1
%
%  combine files from normal COT and moderate SKT  folder


codefilepath = mfilename('fullpath');

% codefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% datafolder, pipelinefolder 
[datafolder, ~, ~, ~] = exp_subfolders();
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


% animal
[fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
animal = codecorresfolder(fi + length('NHPs/'):j);




%%  input setup

% input folder: extracted raw STK data 
inputfolder1 = fullfile(codecorresParentfolder, 'm0_SKTData_extract');
inputfolder2 = fullfile(codecorresParentfolder, 'm0_normalCOTData_extract');

unwanted_DBS = {};

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'avgArea';

%% start here

files1 = dir(fullfile(inputfolder1, '*.mat'));
files2 = dir(fullfile(inputfolder2, '*.mat'));
files = [files1; files2];

nfiles = length(files);
close all;
for fi = 1 : nfiles
    % wait bar
    if(mod(fi, 10) == 0)
        disp(['Avg area lfp data in file ' num2str(fi) '/' num2str(nfiles)]);
    end
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(fi).name;
    load(fullfile(files(fi).folder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    
    
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
    
    
    
    
    % change STN into stn0-1, stn1-2 et.al
    for i = 1: height(T_DBSchnsarea)
        chn = T_DBSchnsarea(i, :).recordingchn;
        barea = cell2mat(T_DBSchnsarea(i, :).brainarea);
        T_DBSchnsarea(i, :).brainarea = {[lower(barea) num2str(chn-1) '-' num2str(chn)]};
        clear chn
    end
    
    % remove unwanted DBS chns
    rows_unwantedDBS = cellfun(@(x) any(strcmp(x, unwanted_DBS)), T_DBSchnsarea.brainarea);
    T_DBSchnsarea(rows_unwantedDBS, :) = [];
    lfpdata_DBS(rows_unwantedDBS, :, :) = [];
    clear rows_unwantedDBS
    
    
    
    % cat GM and DBS
    lfpdata = cat(1, avglfp_notDBS, lfpdata_DBS);
    T_chnsarea = [T_notDBSchnsarea_new; T_DBSchnsarea];
    T_chnsarea.chni = [1: height(T_chnsarea)]';
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx + tmpn - 1) savefilename_addstr ...
        upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];
    
    save(fullfile(savefolder, savefilename), 'lfpdata',  'T_chnsarea', 'fs', 'T_idxevent');
    
    
    
    clear lfpdata fs T_chnsarea T_idxevent 
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

