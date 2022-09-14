function m1_SKTData_avgArea()
%  copy and selected files from Jo\0_dataPrep\SKT\fs500Hz\m1_SKTData_avgArea
%  ref to the files in Jo\0_dataPrep\SKT\fs500Hz\m2_SKTData_SelectTrials
%
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


%%  input setup
inputfolder = fullfile(codecorresParentfolder, '..', 'm1_SKTData_avgArea');
reffolder = fullfile(codecorresParentfolder, '..', 'm2_SKTData_SelectTrials');


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, 'dir')
    rmdir(savecodefolder,'s');
end
copyfile2folder(codefilepath, savecodefolder);


%% Code start here
files = dir(fullfile(reffolder, '*.mat'));

nfiles = length(files);
for filei = 1 : nfiles
    filename = files(filei).name;
    filename = strrep(filename, 'TrialsWMarkers', 'avgAreaSTKData');
    
    success = copyfile(fullfile(inputfolder, filename), fullfile(savefolder, filename));
    
    if ~success
        disp([filename ' not copied']);
    end
end













