function m2_dateBlocks_wUsedM1PMCChns()
%% extract date block pairs for both rest and SKT tasks with used M1 and PMC Channels
%
%   Input:
%       BugMasterDatabase.xlsx: for extracting Rest and SKT tasks
%       Bug_dailyDepth.xlsx: for daily M1 and PMC depth
%   
%   Output:
%       dateBlocks_wUsedM1PMC.csv
%       date_wUsedM1PMC.mat




%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');


% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

% datafolder
[datafolder, ~, ~, ~] = exp_subfolders();


% animal
[fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
animal = codecorresfolder(fi + length('NHPs/'):j);


%% input parameters
xlsxfile_master = fullfile(datafolder, animal, 'BugMasterDatabase.xlsx');
depthFile = fullfile(codecorresParentfolder, 'm1_dailyDepth', 'Bug_dailyDepth.xlsx');


%% global parameters

% different Rest/SKB labels in the master sheet
tasks_Rest = {'Resting', 'Rest'};
tasks_SKB = {'SKB', 'Single'};
tasks = [tasks_Rest tasks_SKB];

strformat_date_master = 'mmddyy';
strformat_date_depth = 'yyyymmdd';
strformat_date_save = 'yyyymmdd';


depth_M1Layer5 = [8 16];
depth_PMCLayer5 = [8 16];


%% save parameters
savefolder = codecorresfolder;
savefilename_dateBlocks = ['dateBlocks_w' num2str(depth_M1Layer5(1)) '-' num2str(depth_M1Layer5(2)) 'UsedM1PMC'];
savefile = fullfile(savefolder, savefilename_dateBlocks);



%% code Start Here

%--- extract t_dateUsedM1PMC_normal/mild/moderate  ---%
disp('extracting t_dateUsedM1PMC_normal/mild/moderate')

depthUsed_M1 = depth_M1Layer5;
depthUsed_PMC = depth_PMCLayer5;

t_depthArea_normal = readtable(depthFile, 'sheet', 'normal');
t_depthArea_mild = readtable(depthFile, 'sheet', 'mild');
t_depthArea_moderate = readtable(depthFile, 'sheet', 'moderate');

t_dateUsedM1PMC_normal = dates_has_usedM1PMC(t_depthArea_normal, depthUsed_M1, depthUsed_PMC); 
t_dateUsedM1PMC_mild = dates_has_usedM1PMC(t_depthArea_mild, depthUsed_M1, depthUsed_PMC); 
t_dateUsedM1PMC_moderate = dates_has_usedM1PMC(t_depthArea_moderate, depthUsed_M1, depthUsed_PMC); 


t_master = readtable(xlsxfile_master);


%----  extract t_dateBlocks ---%

disp('extract t_dateBlocks:')

% normal case
disp('... dealing normal....')
pdcond = "normal";
t_dateUsedM1PMC = t_dateUsedM1PMC_normal;
t_dateBlocks = table();
for di = 1: height(t_dateUsedM1PMC)
    dateofexp = datenum(t_dateUsedM1PMC{di, 1}{1}, strformat_date_depth);
    
    
    % find the tasks records on dateofexp without DBS operation in t_master
    folderName = [animal '_' datestr(dateofexp, strformat_date_master)];
    rowsMask_date = cellfun(@(x) strcmp(x, folderName), t_master.OutputFolderName);
    rowsMask_Task = cellfun(@(x)  any(strcmp(x, tasks)), t_master.BriefDescription);
    rowsMask_noDBS = cellfun(@(x) isempty(x), t_master.DBS_Target);
    rowsMask =  rowsMask_date & rowsMask_Task & rowsMask_noDBS;
    t_tasks_dateofexp = t_master(rowsMask, :);
    
    if isempty(t_tasks_dateofexp)
        continue;
    end
    
    rest_tdtbk = []; rest_mabk = []; 
    SKT_tdtbk = []; SKT_mabk = [];
    for ti = 1: height(t_tasks_dateofexp)
        tdtbk = t_tasks_dateofexp{ti, 'TDTBlock'};
        mabk = t_tasks_dateofexp{ti, 'MAFile'};
        
        if any(strcmp(t_tasks_dateofexp{ti, 'BriefDescription'}, tasks_Rest))
            rest_tdtbk = [rest_tdtbk, tdtbk];
            rest_mabk = [rest_mabk, mabk];
        end
        
        if any(strcmp(t_tasks_dateofexp{ti, 'BriefDescription'}, tasks_SKB))
            SKT_tdtbk = [SKT_tdtbk, tdtbk];
            SKT_mabk = [SKT_mabk, mabk];
        end

        clear tdtbk mabk
    end
    
    if  length(rest_tdtbk) >=2 || length(SKT_tdtbk) >= 2
        continue
    end
    
    if isempty(rest_tdtbk)
        rest_tdtbk = nan;
        rest_mabk = nan;
    end
    if isempty(SKT_tdtbk)
        SKT_tdtbk = nan;
        SKT_mabk = nan;
    end
    
    t_dateblock = table(string(datestr(dateofexp, strformat_date_save)), pdcond,... 
                               rest_tdtbk, rest_mabk, SKT_tdtbk, SKT_mabk, ...
                               t_dateUsedM1PMC{di, 'chnsUsed_M1'}, t_dateUsedM1PMC{di, 'chnsUsed_PMC'}, ...
                               'VariableNames', {'Date', 'pdCond', 'Rest_tdtbk', 'Rest_mabk', 'SKT_tdtbk', 'SKT_mabk', 'chnsUsed_M1', 'chnsUsed_PMC'});
    
                    
    t_dateBlocks = [t_dateBlocks; t_dateblock];
    
    
    
    clear rest_tdtbk rest_mabk SKT_tdtbk SKT_mabk t_dateblock
    clear dateofexp folderName rowsMask_date rowsMask_Task rowsMask_noDBS rowsMask
    clear t_tasks_dateofexp
end
t_dateBlocks_normal =  t_dateBlocks;
clear t_dateBlocks t_dateUsedM1PMC



% mild case
disp('... dealing mild....')
pdcond = "mild";
t_dateUsedM1PMC = t_dateUsedM1PMC_mild;
t_dateBlocks = table();
for di = 1: height(t_dateUsedM1PMC)
    dateofexp = datenum(t_dateUsedM1PMC{di, 1}{1}, strformat_date_depth);
    
    
    % find the tasks records on dateofexp without DBS operation in t_master
    folderName = [animal '_' datestr(dateofexp, strformat_date_master)];
    rowsMask_date = cellfun(@(x) strcmp(x, folderName), t_master.OutputFolderName);
    rowsMask_Task = cellfun(@(x)  any(strcmp(x, tasks)), t_master.BriefDescription);
    rowsMask_noDBS = cellfun(@(x) isempty(x), t_master.DBS_Target);
    rowsMask =  rowsMask_date & rowsMask_Task & rowsMask_noDBS;
    t_tasks_dateofexp = t_master(rowsMask, :);
    
    if isempty(t_tasks_dateofexp)
        continue;
    end
    
    rest_tdtbk = []; rest_mabk = []; 
    SKT_tdtbk = []; SKT_mabk = [];
    for ti = 1: height(t_tasks_dateofexp)
        tdtbk = t_tasks_dateofexp{ti, 'TDTBlock'};
        mabk = t_tasks_dateofexp{ti, 'MAFile'};
        
        if any(strcmp(t_tasks_dateofexp{ti, 'BriefDescription'}, tasks_Rest))
            rest_tdtbk = [rest_tdtbk, tdtbk];
            rest_mabk = [rest_mabk, mabk];
        end
        
        if any(strcmp(t_tasks_dateofexp{ti, 'BriefDescription'}, tasks_SKB))
            SKT_tdtbk = [SKT_tdtbk, tdtbk];
            SKT_mabk = [SKT_mabk, mabk];
        end

        clear tdtbk mabk
    end
    
    if  length(rest_tdtbk) >=2 || length(SKT_tdtbk) >= 2
        continue
    end
    
    if isempty(rest_tdtbk)
        rest_tdtbk = nan;
        rest_mabk = nan;
    end
    if isempty(SKT_tdtbk)
        SKT_tdtbk = nan;
        SKT_mabk = nan;
    end
    
    t_dateblock = table(string(datestr(dateofexp, strformat_date_save)), pdcond,... 
                               rest_tdtbk, rest_mabk, SKT_tdtbk, SKT_mabk, ...
                               t_dateUsedM1PMC{di, 'chnsUsed_M1'}, t_dateUsedM1PMC{di, 'chnsUsed_PMC'}, ...
                               'VariableNames', {'Date', 'pdCond', 'Rest_tdtbk', 'Rest_mabk', 'SKT_tdtbk', 'SKT_mabk', 'chnsUsed_M1', 'chnsUsed_PMC'});
    
                    
    t_dateBlocks = [t_dateBlocks; t_dateblock];
    
    
    
    clear rest_tdtbk rest_mabk SKT_tdtbk SKT_mabk t_dateblock
    clear dateofexp folderName rowsMask_date rowsMask_Task rowsMask_noDBS rowsMask
    clear t_tasks_dateofexp
end
t_dateBlocks_mild =  t_dateBlocks;
clear t_dateBlocks t_dateUsedM1PMC




% moderate case
disp('... dealing moderate....')
t_dateUsedM1PMC = t_dateUsedM1PMC_moderate;
pdcond = "moderate";
t_dateBlocks = table();
for di = 1: height(t_dateUsedM1PMC)
    dateofexp = datenum(t_dateUsedM1PMC{di, 1}{1}, strformat_date_depth);
    
    
    % find the tasks records on dateofexp without DBS operation in t_master
    folderName = [animal '_' datestr(dateofexp, strformat_date_master)];
    rowsMask_date = cellfun(@(x) strcmp(x, folderName), t_master.OutputFolderName);
    rowsMask_Task = cellfun(@(x)  any(strcmp(x, tasks)), t_master.BriefDescription);
    rowsMask_noDBS = cellfun(@(x) isempty(x), t_master.DBS_Target);
    rowsMask =  rowsMask_date & rowsMask_Task & rowsMask_noDBS;
    t_tasks_dateofexp = t_master(rowsMask, :);
    
    if isempty(t_tasks_dateofexp)
        continue;
    end
    
    rest_tdtbk = []; rest_mabk = []; 
    SKT_tdtbk = []; SKT_mabk = [];
    for ti = 1: height(t_tasks_dateofexp)
        tdtbk = t_tasks_dateofexp{ti, 'TDTBlock'};
        mabk = t_tasks_dateofexp{ti, 'MAFile'};
        
        if any(strcmp(t_tasks_dateofexp{ti, 'BriefDescription'}, tasks_Rest))
            rest_tdtbk = [rest_tdtbk, tdtbk];
            rest_mabk = [rest_mabk, mabk];
        end
        
        if any(strcmp(t_tasks_dateofexp{ti, 'BriefDescription'}, tasks_SKB))
            SKT_tdtbk = [SKT_tdtbk, tdtbk];
            SKT_mabk = [SKT_mabk, mabk];
        end

        clear tdtbk mabk
    end
    
    if  length(rest_tdtbk) >=2 || length(SKT_tdtbk) >= 2
        continue
    end
    
    if isempty(rest_tdtbk)
        rest_tdtbk = nan;
        rest_mabk = nan;
    end
    if isempty(SKT_tdtbk)
        SKT_tdtbk = nan;
        SKT_mabk = nan;
    end
    
    t_dateblock = table(string(datestr(dateofexp, strformat_date_save)), pdcond,... 
                               rest_tdtbk, rest_mabk, SKT_tdtbk, SKT_mabk, ...
                               t_dateUsedM1PMC{di, 'chnsUsed_M1'}, t_dateUsedM1PMC{di, 'chnsUsed_PMC'}, ...
                               'VariableNames', {'Date', 'pdCond', 'Rest_tdtbk', 'Rest_mabk', 'SKT_tdtbk', 'SKT_mabk', 'chnsUsed_M1', 'chnsUsed_PMC'});
    
                    
    t_dateBlocks = [t_dateBlocks; t_dateblock];
    
    
    
    clear rest_tdtbk rest_mabk SKT_tdtbk SKT_mabk t_dateblock
    clear dateofexp folderName rowsMask_date rowsMask_Task rowsMask_noDBS rowsMask
    clear t_tasks_dateofexp
end
t_dateBlocks_moderate =  t_dateBlocks;
clear t_dateBlocks t_dateUsedM1PMC



% combine all pdcond
t_dateBlocks = [t_dateBlocks_normal;t_dateBlocks_mild;t_dateBlocks_moderate];




% -- save--  %
save([savefile '.mat'], 't_dateUsedM1PMC_*', 't_dateBlocks')

end % m2_dateBlocks_wUsedM1PMCChns




function t_dateUsedM1PMC = dates_has_usedM1PMC(t_depthArea, depthUsed_M1, depthUsed_PMC) 
%% extract the dates who have both PMC and M1 channels 
%
%   input:
%       t_depthArea: table containing both area and depth, 
%                     output of readtable(depthFile, 'sheet', pdcond)
%       
%       depthUsed_M1, depthUsed_PMC: used depth range for M1 and PMC
%
%
%   output:
%       t_dateUsedM1PMC: table with varNames {'Date', 'chnsUsed_M1', 'chnsUsed_PMC'}

% separate into t_chnArea and t_depth (convert into double)
t_chnArea = t_depthArea(1, :);
depth_array = cellfun(@(x) str2double(x), t_depthArea{2:end, 2:end});
t_depth = [t_depthArea(2:end, 1) array2table(depth_array)];
t_depth.Properties.VariableNames = t_depthArea.Properties.VariableNames;



% t_M1PMCDepth and t_M1PMC_chnArea, containing Date
cols_M1PMC = cellfun(@(x) any(strcmp(x,  {'M1', 'PMC'})), t_chnArea{1,:});
t_M1PMC_Depth = [t_depth(:,1) t_depth(:, cols_M1PMC)];
t_M1PMC_chnArea = [t_chnArea(:, 1) t_chnArea(:,cols_M1PMC)];


% col number of M1 or PMC in t_M1PMC_Depth (same as t_M1PMC_chnArea)
cols_M1 = cellfun(@(x) any(strcmp(x, 'M1')), t_M1PMC_chnArea{1,:});
cols_PMC = cellfun(@(x) any(strcmp(x, 'PMC')), t_M1PMC_chnArea{1,:});
t_dateUsedM1PMC = table();
for di = 1: height(t_M1PMC_Depth) 
    
    t_M1Depth_day = t_M1PMC_Depth(di, cols_M1);
    cols_chnM1Used = find(t_M1Depth_day{:, :} >= depthUsed_M1(1) & t_M1Depth_day{:, :} <= depthUsed_M1(2));
    t_PMCDepth_day = t_M1PMC_Depth(di, cols_PMC);
    cols_chnPMCUsed = find(t_PMCDepth_day{:, :} >= depthUsed_PMC(1) &  t_PMCDepth_day{:, :} <= depthUsed_PMC(2));
    
    if isempty(cols_chnM1Used) || isempty(cols_chnPMCUsed)
        continue;
    end
    
    
    % extract M1 and PMC chn used for each date
    chnsUsed_M1 = []; 
    for ci = 1 : length(cols_chnM1Used)
        chnStr = t_M1Depth_day.Properties.VariableNames{cols_chnM1Used(ci)};
        chn = str2num(chnStr((1+ length('chan') : end)));
        chnsUsed_M1 = [chnsUsed_M1 chn];
        
        clear chnStr chn
    end
    
    chnsUsed_PMC = [];
    for ci = 1 : length(cols_chnPMCUsed)
        chnStr = t_PMCDepth_day.Properties.VariableNames{cols_chnPMCUsed(ci)};
        chn = str2num(chnStr((1+ length('chan') : end)));
        chnsUsed_PMC = [chnsUsed_PMC chn];
        
        clear chnStr chn
    end
    
    t_day = table(t_M1PMC_Depth{di, 1}, {chnsUsed_M1}, {chnsUsed_PMC}, 'VariableNames', {'Date', 'chnsUsed_M1', 'chnsUsed_PMC'});
    t_dateUsedM1PMC = [t_dateUsedM1PMC; t_day];
    
    
    clear t_M1Depth_day cols_chnM1Used t_PMCDepth_day cols_chnPMCUsed 
    clear chnsUsed_M1 chnsUsed_PMC t_day
end
end %dates_has_usedM1PMC



