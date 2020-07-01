function m0_SKTData_extract(animal)
%% extract the STK data  marked by Ying for Pinky
% 
%   1. abandon the file with max trial time > 5s
%
%   2. trial length = max(each trial length) + t_bef + t_aft



if nargin < 1
    animal = 'Bug';
end


%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% codefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% datafolder, pipelinefolder 
[datafolder, ~, ~, ~] = exp_subfolders();
correspipelinefolder = code_corresfolder(codefilepath, true, false);

%% input setup

conds_extract = {'normal', 'mild'};

% channel number of grayMatter
nGM = 96; 

% Input dir:  preprocessed folder in root2
processedfolder_inroot2 = fullfile('/home','lingling','root2','Animals2',animal, 'Recording', 'Processed', 'DataDatabase');

% file data/SKB Beta Analyzed.xlsx for stk information
file_stkInf = fullfile(datafolder, 'SKB Beta Analyzed.xlsx');


% gray matter chn-area information file
filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
file_GMChnsarea =  fullfile(datafolder, filename_GMChnsarea);


%% save setup
savefolder = correspipelinefolder;
savefilename_prefix = [animal '_STKData_'];


%% Starting

tbl_stkInf = readtable(file_stkInf);

% categorical 
tbl_stkInf.Animal = categorical(tbl_stkInf.Animal);
tbl_stkInf.LFPTransfered = categorical(tbl_stkInf.LFPTransfered);
tbl_stkInf.SKB = categorical(tbl_stkInf.SKB);


% extract the rows for Pinky with LFPTransfered == 'Y' & SKB=='Y'
tbl_stkInf = tbl_stkInf(tbl_stkInf.Animal == animal & tbl_stkInf.LFPTransfered == 'Y' & tbl_stkInf.SKB=='Y', :);


% extract all STK trials
close all
f = waitbar(0, ['Extracting all STK trials']);
n = height(tbl_stkInf);
for i = 1 : n
    % waitbar
    waitbar(i/n,f,['i = ' num2str(i) ', Extracting trials in file ' num2str(i) '/' num2str(n)]);
    
    % date of exp, bktdt
    dateofexp = datenum(tbl_stkInf(i,:).Date);
    tdtbk = str2num(char(tbl_stkInf(i,:).TDTBlock));
    
    
    % get the pd conditioon for the date of experiment
    pdcondition = parsePDCondition(dateofexp, animal);
    
    
    % if the pd of the day is not in conds_extract, skip
    if isempty(find(strcmp(pdcondition, conds_extract), 1))
        disp(['pd of ' animal '_' datestr(dateofexp, 'mmddyy') ' is ' pdcondition ', not in conds_extract'])
        
        continue;
    end
    
    % one day path, e.g 
    % processedfolder_inroot2 = '/home/lingling/root2/Animals2/Bug/Recording/Processed/DataDatabase/', 
    % tdtbk = '2', dateofexp = datenum(datetime('2018-06-21'))
    onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);
    
    disp(['extracting ' onedaypath '-tdtbk' num2str(tdtbk)])
    
    
    % extract trials of lfp data for particular date and tdt block number
    [lfptrial_cortical, lfptrial_dbs,fs,T_idxevent, T_dbsChn] = extractlfptrial(onedaypath, tdtbk);
    % skip this day if any is empty 
    if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) ||isempty(T_dbsChn)
        continue;
    end
    

    % concatenate the nGM GM chns and the dbs chns into lfptrial
    lfpdata = cat(1,lfptrial_cortical([1:nGM], :,:),lfptrial_dbs);    
    
    % get the channel inf
    [T_chnsarea] = chanInf(T_dbsChn, nGM, file_GMChnsarea, pdcondition);

    
   
    % save
    savefilename = [savefilename_prefix  pdcondition '_' datestr(dateofexp, 'mmddyy') '_bktdt' num2str(tdtbk)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    clear dateofexp tdtbk onedaypath
    clear lfptrial_cortical lfptrial_dbs fs T_idxevent T_dbsChn chn_cortical
    clear lfpdata T_chnsarea pdcondition savefilename
end

% close the waitbar
close(f);

end






