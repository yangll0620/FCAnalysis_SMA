function sameDates4SKTRestExtract()
%% extract the same dates used for both SKT and Rest

%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% codecorresParentfolder 
[~, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% input setup
animal = 'Pinky';


input_folder_STK = fullfile(codecorresParentfolder, 'SKT', 'm3_STKData_narrowfiltered29_31');
input_folder_Rest = fullfile(codecorresParentfolder, 'Rest', 'm5_restData_segNarrowedDownsample');

%% save setup
savefolder = codecorresParentfolder;
savefilename = 'Pinky_sameDatesUsedforSTKRest';


%% Codes Start Here

files_SKT = dir(fullfile(input_folder_STK, '*.mat'));
a = struct2cell(files_SKT);
filesname_SKT = a(1,:);


files_Rest = dir(fullfile(input_folder_Rest, '*.mat'));
a = struct2cell(files_Rest);
filesname_Rest = a(1,:);


datestrings_rest = {};
datestrings_skt = {};
dates_cond = {};
for fileresti = 1 : length(filesname_Rest)
    
    fileName_rest = filesname_Rest{fileresti};
    
    % find the date string of format yyyymmdd in fileName_rest
    idx = strfind(fileName_rest, '_tdt');
    datestr_rest_yyyymmdd = fileName_rest(idx-8:idx-1);
    
    % convert yyyymmdd to mmddyy used in fileName_skt
    datestr_rest_mmddyy = datestr(datenum(datestr_rest_yyyymmdd, 'yyyymmdd'), 'mmddyy');
    
    % check datestr_skt exist for each cell of filesname_SKT
    if sum(contains(filesname_SKT, datestr_rest_mmddyy)) > 0 % this date exist in both Rest and SKT
        
        datestrings_rest = [datestrings_rest; ['_' datestr_rest_yyyymmdd '_']];
        
        datestrings_skt = [datestrings_skt; ['_' datestr_rest_mmddyy '_']];
        
        
        cond = parsePDCondition(datenum(datestr_rest_mmddyy, 'mmddyy'), 'Pinky');
        dates_cond = [dates_cond; cond];
    end
    
    clear fileName_rest idx datestr_rest_yyyymmdd datestr_rest_mmddyy
end

% create table
tbl_samedates = table(datestrings_rest, datestrings_skt, dates_cond);

% save
writetable(tbl_samedates, fullfile(savefolder, [savefilename '.csv']));









