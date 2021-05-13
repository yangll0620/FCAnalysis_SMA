function m1_restData_cleaned_extract()
% combine all extracted files of different conditions into one

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

%% input setup

inputfolder_normal = fullfile(codecorresParentfolder,'m1_restData_cleaned_extract_normal');
inputfolder_mild = fullfile(codecorresParentfolder,'m1_restData_cleaned_extract_mild3chns');


%%
savefolder = codecorresfolder;

%% code start here
inputfolder = inputfolder_normal;
files = dir(fullfile(inputfolder, '*.mat'));
for fi = 1 : length(files)
    copyfile(fullfile(inputfolder, files(fi).name), fullfile(savefolder, files(fi).name));
end

inputfolder = inputfolder_mild;
files = dir(fullfile(inputfolder, '*.mat'));
for fi = 1 : length(files)
    copyfile(fullfile(inputfolder, files(fi).name), fullfile(savefolder, files(fi).name));
end

