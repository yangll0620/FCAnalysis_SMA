function m0_restDataGM_extract()
%% extract the rest data with same dateofexp in skb analysis in area SMA

%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

[datafolder, ~, pipelinefolder, ~] = exp_subfolders();
correspipelinefolder = code_corresfolder(codefilepath, true, false);

%% input setup
animal = 'Pinky';

% Input dir:  preprocessed folder in root2
processedfolder_inroot2 = ['/home/lingling/root2/Animals2/' animal '/Recording/Processed/DataDatabase'];

% read Pinky_restinf.csv
file_restinf= fullfile(pipelinefolder, ['NHP_' animal], '0_dataPrep', '0_restinf', [animal '_restinf.csv']);
data_restinf = readtable(file_restinf);

% read channel numbers for different area 
filename_chnsarea = fullfile(datafolder, [animal '_GMChnAreaInf.csv']);
T_chnsarea = chnsarea_extract(filename_chnsarea);
chns_recording = T_chnsarea.recordingchn;


%% save setup
savefolder = correspipelinefolder;
savefilename_prefix = [animal '_restdataGM_'];

f = waitbar(0, ['Extracting']);
n =  length(data_restinf.dateofexp);
for i = 1: n
    waitbar(i/n,f,['Extracting ' num2str(i) '/' num2str(n)]);
    
    % date of exp, bktdt, pdcond for one rest
    dateofexp = datetime(num2str(data_restinf.dateofexp(i)), 'InputFormat','yyMMdd');
    bktdt = data_restinf.bktdt(i);
    pdcond = data_restinf.pdCondition{i};
    
    nexfolder_restlfp = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')],...
        'LFP', ['Block-' num2str(bktdt)]);
    
    % extract the prefix of the files and chns in folder_restlfp
    [file_prefix, ~] = nexfilenames_infolder(nexfolder_restlfp);
    
    % extract the lfp data in lSMA and rSMA
    [lfpdata, fs] = lfpdataextract_fromnexfileofallchns(nexfolder_restlfp, file_prefix, chns_recording);
    
    if ~isempty(lfpdata)
        savefile = fullfile(savefolder, [savefilename_prefix  pdcond '_' ...
            datestr(dateofexp, 'mmddyy') '_bktdt' num2str(bktdt)]);
        save(savefile, 'lfpdata','fs', 'T_chnsarea');
    end
 
    clear dateofexp bktdt pdcond nexfolder_restlfp file_prefix chns i
end
close(f);
disp(['extracted lfpdata using grayMatter are saved to ' savefolder])


function T_chnsarea = chnsarea_extract(filename)
%%  read chns_area information from filename
%

T = readtable(filename);
idx = 1;
T_chnsarea = table;
for i = 1 : length(T.brainarea)   
    brainarea = T.brainarea{i};
    
    tmpcell = split(T.channels{i}, ',');
    for j = 1 : length(tmpcell)
        chn = str2num(char(tmpcell{j}));
        newRow = {idx, brainarea, chn};
        T_chnsarea = [T_chnsarea;newRow];
        
        idx = idx + 1;
    end       
end
T_chnsarea.Properties.VariableNames = {'chni', 'brainarea', 'recordingchn'};