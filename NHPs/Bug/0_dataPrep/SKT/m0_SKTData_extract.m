function m0_SKTData_extract()
%% extract the STK data
%
%	Inputs:
%
%		xlsxfile_master
%
%   Steps:
%       1. extract trials from processed folder in server (call function _extractlfptrial
%           a. trial length = max(each trial length) + t_bef + t_aft
%                   t_bef: the time before target on (default: t_bef = 1)
%                   t_aft: the time after mouth (default: t_aft = 0.5)
%           b. only remain the trials markedin both good reach and good return
%
%
%       2. bipolar for DBS channels
%
%       3. Down sample trials into fs_new = 500
%
%       4. add variable T_chnsarea




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


% animal
[i, j] = regexp(correspipelinefolder, 'NHPs/[A-Za-z]*');
animal = correspipelinefolder(i + length('NHPs/'):j);

%% input setup

% channel number of grayMatter
nGM = 96;

% Input dir:  preprocessed folder in root2
processedfolder_inroot2 = fullfile('/home','lingling','root2','Animals2',animal, 'Recording', 'Processed', 'DataDatabase');


% master sheet
xlsxfile_master = fullfile(datafolder, animal, [animal 'MasterDatabase.xlsx']);
strformat_date_master = 'mmddyy'; % the format of date string in master sheet, e.g '012317'


% gray matter chn-area information file
filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
file_GMChnsarea =  fullfile(datafolder, filename_GMChnsarea);




% downsample
fs_new = 500;

% conds
conds_cell = {'normal', 'mild', 'moderate'};


%% save setup
savefolder = correspipelinefolder;
savefilename_prefix = [animal '_STKData_'];


strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'


%% Starting

% ---- extract t_SKT containing the records for STK ---- %
% extract master table
t_master = readtable(xlsxfile_master);


% find the row indices for task_SKB, idx_rows: 0 for not task_SKB, 1 for task_SKB
idx_rows = cellfun(@(x) ~isempty(find(strcmp(x, tasks_SKB))), t_master.BriefDescription);

% table for rows marked SKB labels with only OutputFolderName and TDTBlock columns
t_SKT = t_master(idx_rows, {'OutputFolderName', 'TDTBlock'});

% convert Table Variables from Cell Arrays of Character Vectors to string/double Arrays
t_SKT.OutputFolderName = string(t_SKT.OutputFolderName);
t_SKT.TDTBlock = double(t_SKT.TDTBlock);




% ---- extract all STK trials ---- %
close all
f = waitbar(0, ['Extracting all STK trials']);
nrecords = height(t_SKT);
for i = 1 : nrecords
    % waitbar
    waitbar(i/n,f,['i = ' num2str(i) ', Extracting trials in file ' num2str(i) '/' num2str(n)]);
    
    
    %%% meta data for current record %%%%
    
    % date of exp, bktdt
    outputfoldername = split(t_SKT.OutputFolderName(i), '_');
    dateofexp = datenum(outputfoldername{end}, strformat_date_master);
    tdtbk = t_SKT.TDTBlock(i);
    
    
    % get the pd conditioon for the date of experiment
    pdcondition = parsePDCondition(dateofexp, animal);
    
    
    if ~any(strcmp(pdcondition, conds_cell)) % avoid tomild, tomoderate
        continue;
    end
    
    % one day path
    onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);
    
    disp(['extracting ' onedaypath '-tdtbk' num2str(tdtbk)])
   
    % extract trials of lfp data for particular date and tdt block number
    [lfptrial_cortical, lfptrial_dbs,fs,T_idxevent, T_dbsChn] = extractlfptrial_(onedaypath, tdtbk);
    
    % skip this day if any is empty 
    if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) ||isempty(T_dbsChn)
        continue;
    end
    
    
    
    %%% extract the lfptrial_m1 marked with good channels %%%
    
    
    
    %%% bipolar dbs %%%%
    
    
    %%% down sample %%%
    
    
    
    

    % concatenate the nGM GM chns and the dbs chns into lfptrial
    lfpdata = cat(1,lfptrial_cortical([1:nGM], :,:),lfptrial_dbs);    
    
    % get the channel inf
    [T_chnsarea] = chanInf(T_dbsChn, nGM, file_GMChnsarea, pdcondition);

    
   
    
    
    %%%  save  %%%
    savefilename = [savefilename_prefix  pdcondition '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    
    
    clear dateofexp tdtbk onedaypath
    clear lfptrial_cortical lfptrial_dbs fs T_idxevent T_dbsChn chn_cortical
    clear lfpdata T_chnsarea pdcondition savefilename
end

% close the waitbar
close(f);

end



function [lfptrial_cortical, lfptrial_dbs, fs ,idxeventtbl, chantbl_dbs] = extractlfptrial_(onedaypath, tdtblock)
% extractlfptrial extract trials for LFP data
%
%  [lfptrial_cortical, lfptrial_dbs, chantbl_cortical, chantbl_dbs] =
%  extractlfptrial(onedaypath, block) return extracted LFP trials of 
%  cortical/subcortical data, dbs data, cortical/subcortical channel
%  information and dbs channel information tables, only trials with both
%  only good reach and return are returned
%    
%  Example usage: 
%   onedaypath = '/home/lingling/root2/Animals2/Bug/Recording/Processed/DataDatabase/Bug_062118'
%   tdtblock = 2
%   [lfptrial_cortical, lfptrial_dbs, idxtbl_event,chantbl_cortical, chantbl_dbs] = extractlfptrial(onedaypath, block)      
%
%  Inputs:
%   onedaypath:  'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417'
%   block: 1 
%
%  Used files: Bug\Recording\Processed\DataDatabase\Bug_062118
%       lfpfile_GM      -   \LFP\Block-2\Bug-180621_Block-2_LFPch*.nex
%       lfpfile_dbs     -   \DBSLFP\Block-2\Bug-180621_Block-2_DBSLFP.nex
%       mafile          -   \Block-2\Bug_20180621_2_cleaned_MA_SingleTargetKluver_Analyze2.mat  

%          
%  Outputs:
%   lfptrial_cortical: lfp trials of cortical/subcortical channels 
%                      [chn_cortical * n_temporal * n_trial]
%
%        lfptrial_dbs: lfp trials of dbs channels
%                      [chn_dbs * n_temporal * n_trial], 1-8: STN, 9-16:GP
%        fs: sample rate
%
%
%        idxeventtbl: a table describes the index for events of target onset,
%                      reach onset, touch screen, return and mouth in the trial
%
%         chantbl_dbs:  a table describes each dbs channel information
%
%  More Description:
%       trial length = max(each trial length) + t_bef + t_aft
%       t_bef: the time before target on (default: t_bef = 1)
%       t_aft: the time after mouth (default: t_aft = 0.5)
%       one trailis from 't_target - t_bef'  to 't_mouth + t_after'

%% add NexMatablFiles path
addpath(genpath(fullfile('..', '..', 'toolbox', 'NexMatlabFiles')))

%% MA data
% read the MA data
mafolder = fullfile(onedaypath, ['Block-' num2str(tdtblock)]); %  'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\Block-1'
mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));

% ma file does not exist
if isempty(mafilestruct)
    disp([mafolder 'has no Analyze2.mat'])
    
    lfptrial_cortical = []; lfptrial_dbs = [];
    fs = [];
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

% load SingleTargetKluverMAData
load(fullfile(mafolder, mafilestruct.name)); 

% ma sample rate
fs_ma = SingleTargetKluverMAData.SR;

% time indices for target onset, reach onset, touch screen, return and mouth
TargetTime = SingleTargetKluverMAData.TargetTime;
ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
MouthTimeix = SingleTargetKluverMAData.MouthTimeix;

timeixtbl_ma = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];

% the tag of good reach trials
tag_goodreach = SingleTargetKluverMAData.goodix_reach;
% the tag of good return trials
tag_goodreturn = SingleTargetKluverMAData.goodix_return;

% extract indices of good trials (both have good reach and return)
idx_goodtrials = find(tag_goodreach .* tag_goodreturn == 1);

% if no good trials can be found
if isempty(idx_goodtrials)
    disp(['no good trials are found in ' mafilestruct.name])
    
    lfptrial_cortical = []; lfptrial_dbs = []; 
    fs = []; 
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

%
timeixtbl_ma = timeixtbl_ma(idx_goodtrials,:);

clear TargetTime ReachTimeix TouchTimeix ReturnTimeix MouthTimeix


%% LFP data 
t_bef = 1; t_aft = 0.5; % t_bef: time before target on, t_aft: time after mouth
[n_trial, n_timevars] = size(timeixtbl_ma);

% read each channel data in  LFP data to lfpdata (initialized when is the 1st channel)
folder_cortical = fullfile(onedaypath, 'LFP', ['Block-' num2str(tdtblock)]);

% files in folder_cortical are empty
if isempty(dir(folder_cortical))
    disp([folder_cortical ' has no files!'])
    
    lfptrial_cortical = []; lfptrial_dbs = []; 
    fs = []; 
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

filenames = extractfield(dir(folder_cortical), 'name');
match = cellfun(@(x) ~isempty(regexp(x, ['LFPch[0-9]*.nex'],'match')),filenames,'UniformOutput',false); % match channel file
match = cell2mat(match);
nexnames = filenames(match);


% folder_cortical has no LFPch*.nex files
if isempty(nexnames)
    disp([folder_cortical ' has no .nex files'])
    
        
    lfptrial_cortical = []; lfptrial_dbs = []; 
    fs = []; 
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

% parse and sort the number of channel
for filei = 1 : length(nexnames)
    filename = nexnames{filei};
    tmp = char(regexp(filename, 'ch[0-9]*+.nex', 'match'));
    chns(filei) = str2num(tmp(length('ch')+1: end-length('.nex')));
    
    if filei == 1
        % get the prefix of the each nex file name
        tmp = split(filename, [num2str(chns(filei)) '.nex']);
        file_prefix = tmp{1};
    end
end


chns = sort(chns);
if ~isequal(chns, [1:112])
    disp([folder_cortical 'does not contain 1:112 channels.'])
    
         
    lfptrial_cortical = []; lfptrial_dbs = []; 
    fs = []; 
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

chn_lfp = length(chns);
for i = 1: length(chns)
    filename = [file_prefix num2str(chns(i)) '.nex'];
    [nexlfp_cortical] = readNexFile(fullfile(folder_cortical, filename));
    
    % extract the number of the structure containing LFPchn* data
    name_list = extractfield(cell2mat(nexlfp_cortical.contvars),'name');
    i_lfp = find(contains(name_list, 'LFP')); % % i.e nexlfp_utah.contvars name == 'LFPch1', or 'MUAch1'
    
    if i == 1 % first channel
        fs_lfpcortical = nexlfp_cortical.contvars{i_lfp}.ADFrequency;
        
        % the time index in the LFP neural data based on MA data
        timeixtbl_lfpcortical = timeixtbl_ma;
        timeixtbl_lfpcortical{:,:} = round(timeixtbl_ma{:,:} / fs_ma * fs_lfpcortical);
        
        % initialize lfp_utah : chn_lfp * n_temporal * n_trial
        maxlen = max(timeixtbl_lfpcortical.MouthTimeix - timeixtbl_lfpcortical.TargetTime) + 1; % maximum length across all trials (unit: ind)
        n_bef = round(t_bef * fs_lfpcortical);  % n_bef: index number before target on
        n_aft = round(t_aft * fs_lfpcortical);  % n_aft: index number after mouth
        n_temporal = maxlen + n_bef + n_aft;
        idx_str  = timeixtbl_lfpcortical.TargetTime - n_bef;
        idx_end  = timeixtbl_lfpcortical.TargetTime + (maxlen -1) + n_aft;
        lfptrial_cortical = zeros(chn_lfp, n_temporal, n_trial);
        
        % idxtbl_lfptrialutah: the idx for events of target onset, reach onset, touch screen,
        % return and mouth in the trial matrix lfpdata_utah (chn_lfputah * n_temporal * n_trial)
        % the first sample corresponds to target onset - t_bef
        idxtbl_lfptrial_cortical = timeixtbl_lfpcortical;
        idxtbl_lfptrial_cortical{:,:} = idxtbl_lfptrial_cortical{:,:} - repmat(idx_str,[1, n_timevars]);
        
        clear timeixtbl_lfpsepchn n_bef n_aft 
    else
        if fs_lfpcortical ~= nexlfp_cortical.contvars{i_lfp}.ADFrequency % samping frequency is different
            chni = chns(i);
            disp(['sampling frequency is different for chni = ' num2str(chni)]);
            break;
        end
    end
    
    % extract each trial for lfp data stored in separate channel
    for triali = 1: n_trial
        if size(nexlfp_cortical.contvars{i_lfp}.data,1) < idx_end(triali)
            disp(mafilestruct)
            disp(filename)
            disp(fs_ma)
            disp(fs_lfpcortical)
            disp(['idx in ma' num2str(timeixtbl_ma{1,end})])
            disp(size(nexlfp_cortical.contvars{i_lfp}.data,1))
            disp(idx_end(triali))
            disp(['maxlen  = ' num2str(maxlen)])
        end
        lfptrial_cortical(i, :, triali) = nexlfp_cortical.contvars{i_lfp}.data(idx_str(triali):idx_end(triali));
        
    end
    
    clear filename i_lfp 
end

% disp play the max trial time
disp(['max trial time is ' num2str(maxlen/fs_lfpcortical)]);

if maxlen/fs_lfpcortical > 5
    disp('Abandon the file with max trial time > 5s')
    lfptrial_cortical = []; lfptrial_dbs = []; fs =[]; idxeventtbl = []; chantbl_dbs = [];
    return
end

%% DBSLFP data
% read DBSLFP in nex file
dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtblock)]); % 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\DBSLFP\Block-1'

% dbs file does not exist
dbsfilepattern = fullfile(dbslfpfolder, ['*Block-' num2str(tdtblock) '_DBSLFP.nex']); 
if  isempty(dir(dbsfilepattern))
    disp([dbsfilepattern, ' files do not exit!'])
    
    lfptrial_cortical = []; lfptrial_dbs = [];
    fs = [];
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

filenames = extractfield(dir(dbsfilepattern),'name'); % Bug-180621_Block-2_DBSLFP.nex
match = cellfun(@(x) ~isempty(regexp(x, ['[A-Z, a-z, 0-9]*.nex'],'match')),filenames,'UniformOutput',false); % match channel file
match = cell2mat(match);
nexnames = filenames(match);
if length(nexnames) ~= 1
    % only has one file inside folder
    error(['More than one _DBSLFP.nex  in file' dbslfpfolder])
end
[nexlfp_dbs] = readNexFile(fullfile(dbslfpfolder, nexnames{1}));

% find the STN and GP data in the struct nexlfp_dbs
% extract all the values of field 'name' in nexlfp_dbs.convars
convars = cell2mat(nexlfp_dbs.contvars);
name_list = extractfield(cell2mat(nexlfp_dbs.contvars),'name');
% extract the indices representing STN and GP data in nexlfp_dbs.contvars (33*1 struct array)
idx_dbs = find(contains(name_list, 'RAW_DBSch'));


% check the dbs channel frequencies
if range(cell2mat({convars(idx_dbs).ADFrequency})) ~= 0
    % check the sampling frequencies in channels are consistent
    disp(['ADFrequency in the dbs channels (STN & GP) is not consistent']);
    return
end
if(convars(idx_dbs(1)).ADFrequency ~= fs_lfpcortical) 
    % check the sampling frequency of dbs is the same as the lfp stored in separate channels or not
    disp(['the sampling frequency of dbs is not the same as the lfp stored in separate channels'])
    return
end
fs = fs_lfpcortical; % sampling frequency for dbs channels
idxeventtbl = idxtbl_lfptrial_cortical;
clear idxtbl_lfptrial_sepchn

% extract each trial for lfp dbs data
nchn_dbs = length(idx_dbs);
lfptrial_dbs = zeros(nchn_dbs, n_temporal, n_trial);
for i = 1: nchn_dbs
    chni = idx_dbs(i);
    lfp_1chn = convars(chni).data;
    for triali = 1: n_trial
        lfptrial_dbs(i, :, triali) = lfp_1chn(idx_str(triali):idx_end(triali));
    end
    clear chni lfp_1chn triali
end

% chantbl_dbs: a table describes each dbs channel information
elecchn = extractfield(convars(idx_dbs),'name');
elecchn = elecchn';
area = cell(nchn_dbs, 1);
area(1:8) = {'STN'};
area(9: 16) = {'GP'};
chantbl_dbs =[table(area) table(elecchn)];
clear varName

clear convars filename idx_stn idx_gp 
end






