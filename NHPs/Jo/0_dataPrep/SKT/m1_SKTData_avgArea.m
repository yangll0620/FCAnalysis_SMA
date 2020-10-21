function m1_SKTData_avgArea()
%% average lfp across each GM area
%
%   


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
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');



%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'avgArea';

%% start here

files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
close all;
for fi = 1 : nfiles
    % wait bar
    if(mod(fi, 10) == 0)
        disp(['Avg area lfp data in file ' num2str(fi) '/' num2str(nfiles)]);
    end
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(fi).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    
    
    % average across each GM area
    mask_GM = strcmp(T_chnsarea.electype, 'Gray Matter');
    mask_DBS = strcmp(T_chnsarea.electype, 'DBS');
    lfpdata_GM = lfpdata(mask_GM, :, :);
    lfpdata_DBS = lfpdata(mask_DBS, :, :);
    T_GMchnsarea = T_chnsarea(mask_GM, :);
    T_DBSchnsarea = T_chnsarea(mask_DBS, :);
    [avglfp_GM, T_GMchnsarea_new] = avglfp_acrossGMArea(lfpdata_GM, T_GMchnsarea);
    
    if isempty(avglfp_GM)
        continue;
    end
    
    % cat GM and DBS
    lfpdata = cat(1, avglfp_GM, lfpdata_DBS);
    T_chnsarea = [T_GMchnsarea_new; T_DBSchnsarea];
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

function [avglfp, T_GMchnsarea_new] = avglfp_acrossGMArea(lfpdata_GM, T_GMchnsarea)
% average lfp across each GM area 
% 
% args:
%       lfpdata_GM : lfp data in GM (nchns * ntemp * ntrials )
%       T_GMchnsarea: chns-area table in GM
% return:
%       avglfp: averaged lfp data (nareas * ntemp * ntrials
%       T_GMchnsarea_new: the new chns-area table for GM



uniqGMAreas = unique(T_GMchnsarea.brainarea);
avglfp = [];
T_GMchnsarea_new = T_GMchnsarea([], :);
for ai = 1 : length(uniqGMAreas)
    GMarea = uniqGMAreas{ai};
   
    mask_area = strcmp(T_GMchnsarea.brainarea, GMarea);
    lfp_area = lfpdata_GM(mask_area, :, :);
    avglfp_area = mean(lfp_area, 1);
    
    
    avglfp = cat(1, avglfp, avglfp_area);
    
    T_GMchnsarea_new = [T_GMchnsarea_new; {ai, GMarea, nan, 'Gray Matter', []}];
    
    clear GMarea mask_area avglfp_area
end



end

