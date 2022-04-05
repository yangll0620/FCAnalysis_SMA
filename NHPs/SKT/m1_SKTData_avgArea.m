function m1_SKTData_avgArea(animal)
%% average lfp across each area
%
%
%   1. remove chns marked with multiple areas or GP/STN in Gray Matter
%   2. average across each GM area (no DBS)
%   3. change STN/GP into stn0-1, stn1-2... gp0-1, gp1-2 et.al
%   4. remove unwanted DBS chns



%% params
if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end

codefilepath = mfilename('fullpath');

% codefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path and NHP path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));


[~, codefilename]= fileparts(codefilepath);
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, 'SKT', codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);


%%  input setup

% input folder: extracted raw STK data 
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');


%% save setup
savefolder = codecorresfolder;

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

savefilename_prefix = 'avgArea';

%% start here
unwanted_DBS = unwanted_DBS_extract(animal);
files = dir(fullfile(inputfolder, '*.mat'));


nfiles = length(files);
for fi = 1 : nfiles
    % wait bar
    disp(['Avg area lfp data in file ' num2str(fi) '/' num2str(nfiles)]);
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(fi).name;
    load(fullfile(files(fi).folder, filename), 'lfpdata', 'fs_lfp', 'mask_goodreach', 'mask_goodreturn', 'T_chnsarea', 'T_idxevent_lfp', ...
                                               'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    
    % remove chns marked with multiple areas or GP/STN in Gray Matter
    mask_multipleAreas = cellfun(@(x) contains(x, '/'), T_chnsarea.brainarea);
    mask_GPSTNinGM = cellfun(@(x) contains(x, 'GP'), T_chnsarea.brainarea) & ~strcmp(T_chnsarea.electype, 'DBS');
    T_chnsarea = T_chnsarea(~mask_multipleAreas & ~mask_GPSTNinGM, :); 
    T_chnsarea.chni = [1:height(T_chnsarea)]';
    for tri = 1 : length(lfpdata)
        lfp_1trial = lfpdata{tri};
        lfp_1trial = lfp_1trial(~mask_multipleAreas & ~mask_GPSTNinGM, :);
        lfpdata{tri} = lfp_1trial;
        clear lfp_1trial
    end
    
       
    % average across each notDBS area and create new T_chnsarea
    mask_notDBS = ~strcmp(T_chnsarea.electype, 'DBS');
    uniqAreas = unique(T_chnsarea.brainarea(mask_notDBS));
    T_chnsarea_notDBS = T_chnsarea([], :);
    for tri = 1: length(lfpdata)
        lfp_1trial = lfpdata{tri};
        avglfp_notDBS = [];
        for ai = 1 : length(uniqAreas)
            areaNotDBS = uniqAreas{ai};
            mask_areaNotDBS = mask_notDBS & strcmp(T_chnsarea.brainarea, areaNotDBS);
            lfp_1trial_notDBSArea = lfp_1trial(mask_areaNotDBS, :);
            avglfp_area = mean(lfp_1trial_notDBSArea, 1);
            avglfp_notDBS = cat(1, avglfp_notDBS, avglfp_area);
            
            if tri == 1 % create new T_chnsarea_notDBS
                T_chnsarea_notDBS = [T_chnsarea_notDBS; {ai, areaNotDBS, nan,  'Gray Matter', []}];
            end
            
            clear areaNotDBS mask_areaNotDBS lfp_1trial_notDBSArea avglfp_area
        end
        
        lfpdata{tri} = cat(1, avglfp_notDBS, lfp_1trial(~mask_notDBS, :));
        clear lfp_1trial avglfp
    end
    T_chnsarea = [T_chnsarea_notDBS; T_chnsarea(~mask_notDBS, :);];
    T_chnsarea.chni = [1: height(T_chnsarea)]';
    
    
    
    % change STN/GP into stn0-1, stn1-2... gp0-1, gp1-2 et.al
    for i = 1: height(T_chnsarea)
        barea = cell2mat(T_chnsarea(i, :).brainarea);
        if ~(strcmpi(barea, 'STN') || strcmpi(barea, 'GP'))
            continue;
        end
        chn = T_chnsarea(i, :).recordingchn;
        T_chnsarea(i, :).brainarea = {[lower(barea) num2str(chn-1) '-' num2str(chn)]};
        clear chn
    end
    
    
    % remove unwanted DBS chns
    rows_unwantedDBS = cellfun(@(x) any(strcmp(x, unwanted_DBS)), T_chnsarea.brainarea) & strcmp(T_chnsarea.electype, 'DBS');
    T_chnsarea(rows_unwantedDBS, :) = [];
    T_chnsarea.chni = [1: height(T_chnsarea)]';
    for tri = 1: length(lfpdata)
        lfp_1trial = lfpdata{tri};
        lfp_1trial(rows_unwantedDBS, :) = [];
        lfpdata{tri} = lfp_1trial;
        clear lfp_1trial
    end
    clear rows_unwantedDBS
    
    
    % save
    idx = strfind(filename, [animal]);
    tmpn = length([animal]);
    savefilename = [filename(idx:idx + tmpn - 1) '_' savefilename_prefix ...
        upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];
    clear idx tmpn
    
    save(fullfile(savefolder, savefilename), 'lfpdata', 'fs_lfp', 'mask_goodreach', 'mask_goodreturn','T_chnsarea', 'T_idxevent_lfp', ...
        'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    
    
    
    clear lfpdata fs_lfp mask_goodReach mask_goodReturn T_chnsarea T_idxevent_lfp fs_ma  T_idxevent_ma smoothWspeed_trial Wpos_smooth_trial Wrist_smooth_trial
    clear savefilename
end



end
