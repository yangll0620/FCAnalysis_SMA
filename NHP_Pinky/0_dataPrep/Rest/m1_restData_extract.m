function m1_restData_extract()
%% extract the rest data with same dateofexp in skb analysis in all recorded Channel


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
[datafolder, ~, pipelinefolder, ~] = exp_subfolders();
[corrPPLfolder, corrPPLParfolder] = code_corresfolder(codefilepath, true, false);

%% input setup
animal = 'Pinky';

% Input dir:  preprocessed folder in root2
processedfolder_inroot2 = ['/home/lingling/root2/Animals2/' animal '/Recording/Processed/DataDatabase'];

% STK trial Information file name
filefolder_restInf = fullfile(corrPPLParfolder, 'm0_restinf');
file_restInf = fullfile(filefolder_restInf, [animal '_restinf.csv']);


% GrayMatter chn-area information file
filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
file_GMChnsarea =  fullfile(datafolder, filename_GMChnsarea);


%% save setup
savefolder = corrPPLfolder;
savefilename_prefix = [animal '_restData_'];



%% Starting
% Read stk information into table
metadata_restInf = readtable(file_restInf);

% extract all STK trials
f = waitbar(0, ['Extracting all rest Data with same dateofexp in skb analysis']);
nRows = height(metadata_restInf);

% channel number of utah array, grayMatter, and dbs leads
nUtah = 96; nGM = 32; nDBS = 16;
f = waitbar(0, ['Extracting rest Data....']);
for rowi = 1: nRows
    
    % waitbar
    waitbar(rowi/nRows,f,['Extracting rest Data in file ' num2str(rowi) '/' num2str(nRows)]);
    
    % date of exp, bktdt, pdcond for one rest
    dateofexp = datetime(num2str(metadata_restInf.dateofexp(rowi)), 'InputFormat','yyMMdd');
    bktdt = metadata_restInf.bktdt(rowi);
    pdcond = metadata_restInf.pdCondition{rowi};
    
    
    % one day path
    onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);
    
    disp(['extracting ' onedaypath '-tdtbk' num2str(bktdt)])
 
    
     % extract rest of lfp data for particular date and tdt block number
    [lfpRest_UtahGM, lfpRest_DBS, fs, tbl_ChanDBS] = extractRestLFP(onedaypath, bktdt);
    
    if ~isempty(lfpRest_UtahGM)
        
        % get the channel inf
        [T_chnsarea] = chanInf(tbl_ChanDBS, nUtah, nGM, file_GMChnsarea, nDBS);
        
        % concatenate the lfp data of utah, GM and dbs chns into lfptrial
        lfpdata = cat(1,lfpRest_UtahGM, lfpRest_DBS);
        
        % save
        savefile = fullfile(savefolder, [savefilename_prefix  pdcond '_' ...
            datestr(dateofexp, 'mmddyy') '_bktdt' num2str(bktdt)]);
        save(savefile, 'lfpdata', 'fs', 'T_chnsarea');
        
        clear tbl_ChnsArea lfpdata savefile
    end
    
    
    clear dateofexp bktdt pdcond onedaypath 
    clear lfpRest_UtahGM lfpRest_DBS fs tbl_ChanDBS
end
close(f);
disp(['Extraced rest Data in ' savefolder])
end



function [tbl_ChnsArea] = chanInf(tbl_ChanDBS, nUtah, nGM, file_GMChnsarea ,nDBS)
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
brainareas(nUtah+nGM+1:nUtah+nGM+nDBS,1) = tbl_ChanDBS.area;
notes(nUtah+nGM+1:nUtah+nGM+nDBS,1) = tbl_ChanDBS.elecchn;
recordingchn(nUtah+nGM+1:nUtah+nGM+nDBS) = [1:nDBS];

% channel information table
tbl_ChnsArea = table;
tbl_ChnsArea.chni = chni_vec';
tbl_ChnsArea.brainarea = brainareas;
tbl_ChnsArea.recordingchn = recordingchn;
tbl_ChnsArea.electype = electype;
tbl_ChnsArea.notes = notes;
end