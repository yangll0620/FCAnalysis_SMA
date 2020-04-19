function m0_STKData_extract()
%% extract the STK data  marked by Ying for Pinky


%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% datafolder, pipelinefolder 
[datafolder, ~, ~, ~] = exp_subfolders();
correspipelinefolder = code_corresfolder(codefilepath, true, false);

%% input setup
animal = 'Pinky';

% Input dir:  preprocessed folder in root2
processedfolder_inroot2 = ['/home/lingling/root2/Animals2/' animal '/Recording/Processed/DataDatabase'];

% STK trial Information file name
filename_stkInf = [animal '_STK_markedbyYing.csv'];
file_stkInf = fullfile(datafolder, filename_stkInf);


% gray matter chn-area information file
filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
file_GMChnsarea =  fullfile(datafolder, filename_GMChnsarea);


%% save setup
savefolder = correspipelinefolder;
savefilename_prefix = [animal '_STKData_'];


%% Starting

% Read stk information into table
metadata_stkInf = readtable(file_stkInf);
% extract the dateofexp, bkma, and bktdt of used STK trials which are marked 'Yes' in the column 'YingUsed'
validstks = metadata_stkInf{strcmp(metadata_stkInf.YingUsed, 'Yes'), {'dateofexp','bkma', 'bktdt'}};


% extract all STK trials
f = waitbar(0, ['Extracting all STK trials']);
n = length(validstks);
% channel number of utah array, grayMatter, and dbs leads
nUtah = 96; nGM = 32; nDBS = 16;
for i = 1 : n
    % waitbar
    waitbar(i/n,f,['Extracting trials in file ' num2str(i) '/' num2str(n)]);
    
    % date of exp, bktdt
    dateofexp = datenum(string(validstks(i, 1)), 'yymmdd');
    bktdt = validstks(i, 3);
    
    % one day path
    onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);
    
    disp(['extracting ' onedaypath '-tdtbk' num2str(bktdt)])
    
    % extract trials of lfp data for particular date and tdt block number
    [lfptrial_cortical, lfptrial_dbs,fs,T_idxevent, T_dbsChn] = extractlfptrial(onedaypath, bktdt);
    
    % skip this day if any is empty 
    if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) ||isempty(T_dbsChn)
        continue;
    end
    
    
    chn_cortical = size(lfptrial_cortical, 1);
    % concatenate the 96 chns utah, 32 chns GM and the dbs chns into lfptrial
    lfpdata = cat(1,lfptrial_cortical([1:96 chn_cortical-31:chn_cortical], :,:),lfptrial_dbs);
    
    % get the channel inf
    [T_chnsarea] = chanInf(T_dbsChn, nUtah, nGM, file_GMChnsarea, nDBS);
    
    % get the pd conditioon for the date of experiment
    pdcondition = parsePDCondition_Pinky(dateofexp);
    
    
    % save
    savefilename = [savefilename_prefix  pdcondition '_' datestr(dateofexp, 'mmddyy') '_bktdt' num2str(bktdt)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    clear dateofexp bktdt onedaypath
    clear lfptrial_cortical lfptrial_dbs fs T_idxevent T_dbsChn chn_cortical
    clear lfpdata T_chnsarea pdcondition savefilename
end
end


function [T_chnsarea] = chanInf(chantbl_dbs, nUtah, nGM, file_GMChnsarea ,nDBS)
% extract channel inf table
%
%   Output:
%       T_chnsarea: table of channel inf

chni_vec = uint8([1:nUtah+nGM+nDBS]);
electype = cell(nUtah+nGM+nDBS,1);
brainareas = cell(nUtah+nGM+nDBS,1);
notes = cell(nUtah+nGM+nDBS,1);
recordingchn = zeros(nUtah+nGM+nDBS,1);

% electype
electype(1:nUtah,1) = {'Utah Array'};
electype(nUtah+1:nUtah+nGM,1) = {'Gray Matter'};
electype(nUtah+nGM+1:nUtah+nGM+nDBS,1) = {'DBS'};

% deal with GM, all chns of utah array area for M1 
brainareas(1:nUtah,1) = {'M1'};
notes(1:nUtah,1) = {''};
recordingchn(1:nUtah) = [1:nUtah];

% deal with Gray Matter
brainareas(nUtah+1:nUtah+nGM,1) = {''};
T = readtable(file_GMChnsarea);
chi_firstGM = 101;
for i = 1 : length(T.brainarea)
    area = T.brainarea{i};
    tmpcell = split(T.channels{i}, ',');
    
    for j = 1 : length(tmpcell)
        chn = str2num(char(tmpcell{j}));
        brainareas(chn-chi_firstGM+nUtah+1,1) = {area};
    end
end
recordingchn(nUtah+1:nUtah+nGM) = [nUtah+1:nUtah+nGM] + (chi_firstGM - nUtah-1);
notes(nUtah+1:nUtah+nGM,1) = {''};

% deal with the DBS channel table
brainareas(nUtah+nGM+1:nUtah+nGM+nDBS,1) = chantbl_dbs.area;
notes(nUtah+nGM+1:nUtah+nGM+nDBS,1) = chantbl_dbs.elecchn;
recordingchn(nUtah+nGM+1:nUtah+nGM+nDBS) = [1:nDBS];

% channel information table
T_chnsarea = table;
T_chnsarea.chni = chni_vec';
T_chnsarea.brainarea = brainareas;
T_chnsarea.recordingchn = recordingchn;
T_chnsarea.electype = electype;
T_chnsarea.notes = notes;
end

