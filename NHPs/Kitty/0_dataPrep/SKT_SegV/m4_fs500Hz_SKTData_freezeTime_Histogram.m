function m4_fs500Hz_SKTData_freezeTime_Histogram()
% Objective:
%       extract freezing episode
%       data structure
%           episodes{}
% 1. manupulation time > 5s

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

inputfolder = fullfile(codecorresParentfolder, 'm3_fs500Hz_SKTData_freezeEpisodeExtract');


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code Start Here
optFreezeTypes = optFreezeTypes_extract('codesavefolder', savecodefolder);
tdur_freeze = zeros(length(optFreezeTypes), 1);

files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
for fi = 1 : length(files)
    filename = files(fi).name;
    
    load(fullfile(inputfolder, filename), 'freezStruct');
    freezEpisodes = freezStruct.freezEpisodes;
    for freezi = 1 : length(freezEpisodes)
        freezeType = freezEpisodes{freezi}.freezeType;
        mask_freeze = strcmp(freezeType, optFreezeTypes);
        tdur_freeze(mask_freeze) = tdur_freeze(mask_freeze) + 
    end
    
    clear filename freezStruct
end 
