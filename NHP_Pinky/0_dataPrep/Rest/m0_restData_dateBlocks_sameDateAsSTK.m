function m0_restData_dateBlocks_sameDateAsSTK()
%% extract the dateBlock pairs for rest data whose dates are the same as the STK data
%
%
%
%   Inputs:
%   
%       sktFolder:   NHP_Pinky/0_dataPrep/SKT/m3_STKData_narrowfiltered29_31 for extracting STK dates
%
%       masterFile:   datafolder/PinkyMasterDatabase.xlsx for extracting corresponding rest tdtblock
%
%   Output:
%       
%       Pinky_restDateBlocks_sameDatesAsSTK.mat: variable restDateBlocks_sameDatesAsSTK

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
[codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

% datafolder
[datafolder, ~, pipelinefolder, ~] = exp_subfolders();

%% input setup

% stk folder for extracting the dates used for STK analysis
sktFolder = fullfile(pipelinefolder, 'NHP_Pinky', '0_dataPrep', 'SKT', 'm3_STKData_narrowfiltered29_31');
stkfiles = dir(fullfile(sktFolder, '*.mat'));

% master file for extracting the corresponding rest blocki for each used date
masterFile = fullfile(datafolder, 'PinkyMasterDatabase.xlsx');



%% save setup
savefolder = codecorresfolder;
savefilename = 'Pinky_restDateBlocks_sameDatesAsSTK.mat';


%% Start Here
tbl_master = readtable(masterFile);

% categorical 
tbl_master.OutputFolderName = categorical(tbl_master.OutputFolderName);
tbl_master.BriefDescription = categorical(tbl_master.BriefDescription);

restDateBlocks_sameDatesAsSTK = {};
for i = 1: length(stkfiles)
    filename = getfield(stkfiles, {i}, 'name');
    
    idx = strfind(filename, '_bktdt');
    datestring = filename(idx-6 : idx-1);
    
    
    % extract the rows for Pinky rest in date using column OutputFolderName and BriefDescription
    tbl_date = tbl_master(tbl_master.OutputFolderName == ['Pinky_' datestring]  & tbl_master.BriefDescription=='Resting', :);
    
    % combine the extract resting dateblock string into restDateBlocks_sameDatesAsSTK
    for rowi = 1 : height(tbl_date)
        tdtblock = getfield(tbl_date, {rowi}, 'TDTBlock');
        
        % generate the dateBlockString e.g. '20170915_1'
        dateblockstring = [datestr(datenum(datestring, 'mmddyy'), 'yyyymmdd') '_' num2str(tdtblock)];
        
        % combine the new obtained datestr
        restDateBlocks_sameDatesAsSTK = [restDateBlocks_sameDatesAsSTK; dateblockstring];
        
        clear tdtblock dateblockstring
    end
    

    
    
    clear filename idx datestr tbl_date sktUsedDates rowi
end

save(fullfile(savefolder, savefilename), 'restDateBlocks_sameDatesAsSTK');