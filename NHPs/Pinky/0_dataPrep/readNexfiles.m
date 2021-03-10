nexfilefolder = '/home/lingling/root2/Animals2/Pinky/Recording/Processed/DataDatabase/Pinky_103017/LFP/Block-1';

chns_M1 = [1:96];

%% add util folder to path
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));


%% exp_subfolders & correfolder
[datafolder, codefolder, pipelinefolder, outputfolder] = exp_subfolders();

codecorresfolder = code_corresfolder(codefilepath);

%% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))




%% read
[datnex] = readNexFile(nexfile);