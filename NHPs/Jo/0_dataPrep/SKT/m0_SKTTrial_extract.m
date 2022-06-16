function m0_SKTTrial_extract()
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
%       3. add variable T_chnsarea


%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1:strfind(codefilepath, 'code') + length('code') - 1);

% add util path
addpath(genpath(fullfile(codefolder, 'util')));
addpath(genpath(fullfile(codefolder, 'NHPs')));
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))



% datafolder, pipelinefolder
[datafolder, ~, ~, ~] = exp_subfolders();
correspipelinefolder = code_corresfolder(codefilepath, true, false);

%% input setup
animal = animal_extract(correspipelinefolder);

% Input dir:  preprocessed folder in root2
processedfolder_inroot2 = fullfile('H:', 'My Drive', 'NMRC_umn', 'Projects', 'FCAnalysis', 'exp',  'data', 'Animals', animal, 'Recording', 'Processed', 'DataDatabase');

% master sheet
xlsxfile_master = fullfile(datafolder, 'Animals', animal, [animal 'MasterDatabase.xlsx']);
strformat_date_master = 'mmddyy'; % the format of date string in master sheet, e.g '012317'

% different SKB labels in the master sheet
tasks_SKB = {'SKB', 'Single Target Kluver', 'Single target Kluver', 'Single target kluver', 'Stingle target Kluver', 'single target kluver'};


% channel number of utah array, grayMatter, and dbs leads
nM1 = 96; nSTN = 8; nGP = 8;


% conds
conds_cell = cond_cell_extract(animal);

%% save setup
savefolder = correspipelinefolder;
savecodefolder = fullfile(correspipelinefolder, 'code');
savefilename_prefix = [animal '_STKData_'];

strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'


%% Starting Here

if ~exist(savecodefolder, 'dir')
    mkdir(savecodefolder);
end
copyfile2folder(codefilepath, savecodefolder);

% extract master table
t_master = readtable(xlsxfile_master);

% find the row indices for task_SKB, idx_rows: 0 for not task_SKB, 1 for task_SKB
idx_rows = cellfun(@(x) ~isempty(find(strcmp(x, tasks_SKB))), t_master.BriefDescription);

% table for rows marked SKB labels with only OutputFolderName and TDTBlock columns
t_SKT = t_master(idx_rows, {'OutputFolderName', 'TDTBlock'});

% convert Table Variables from Cell Arrays of Character Vectors to string/double Arrays
t_SKT.OutputFolderName = string(t_SKT.OutputFolderName);
t_SKT.TDTBlock = double(t_SKT.TDTBlock);

%%% add chn-area information T_chnsarea  %%%
T_chnsarea = chanInf_M1DBS(nM1, nSTN -1, nGP -1); % bipolar DBS


%%% extract all STK trials %%%
f = waitbar(0, 'Extracting all STK trials');
nrecords = height(t_SKT);

for i = 1:nrecords
    % waitbar
    waitbar(i / nrecords, f, ['Extracting trials in file ' num2str(i) '/' num2str(nrecords)]);
    
    % date of exp, bktdt
    outputfoldername = split(t_SKT.OutputFolderName(i), '_');
    dateofexp = datenum(outputfoldername{end}, strformat_date_master);
    tdtbk = t_SKT.TDTBlock(i);
    
    
    % get the pd conditioon for the date of experiment
    pdcond = parsePDCondition(dateofexp, animal);
    
    
    if ~any(strcmp(pdcond, conds_cell)) % avoid tomild, tomoderate
        
        clear outputfoldername dateofexp tdtbk
        clear pdcond
        
        continue;
    end
    
    
    savefilename = [savefilename_prefix pdcond '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk) '.mat'];
    if exist(fullfile(savefolder, savefilename), 'file') == 2
        clear outputfoldername dateofexp tdtbk
        clear pdcond savefilename
        
        continue;
    end
    
    % one day path
    onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);
    if ~exist(onedaypath, 'dir')
        clear outputfoldername dateofexp tdtbk
        clear pdcond savefilename onedaypath
        continue;
    end
    
    disp(['extracting ' onedaypath '-tdtbk' num2str(tdtbk)])
    
    % extract trials of lfp data for particular date and tdt block number
    [lfptrial_notDBS, lfptrial_dbs, fs_lfp, T_idxevent_lfp, fs_ma, T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial, T_dbsChn] = ...
        extract_lfptrial(onedaypath, tdtbk);
    
    % skip this day if any is empty
    if isempty(lfptrial_notDBS) || isempty(lfptrial_dbs) || isempty(fs_lfp) || isempty(T_idxevent_lfp) || isempty(T_dbsChn)
        clear outputfoldername dateofexp tdtbk
        clear pdcond savefilename onedaypath
        clear lfptrial_cortical lfptrial_dbs fs_lfp T_idxevent_lfp
        clear fs_ma T_idxevent_ma smoothWspeed_trial Wpos_smooth_trial Wrist_smooth_trial T_dbsChn
        
        continue;
    end
    
    
    %%% extract the lfptrial_m1 marked with good channels %%%
    if strcmpi(animal, 'Jo')
        lfptrial_m1 = lfptrial_notDBS;
    end
    
    %%% bipolar dbs %%%%
    for tri = 1 : length(lfptrial_dbs)
        tmp = lfptrial_dbs{tri}; % tmp: 16 * ntemp
        bipo_stn = diff(tmp(1:nSTN, :), [], 1);
        bipo_gp = diff(tmp(nSTN +1:nSTN+nGP, :), [], 1);
        lfptrial_dbs{tri} = cat(1, bipo_stn, bipo_gp);
        clear tmp bipo_stn bipo_gp
    end

    
    
    % concatenate the utah chns, and the dbs chns into lfptrial
    for tri = 1 : length(lfptrial_dbs)
        lfpdata{tri} = cat(1, lfptrial_m1{tri}, lfptrial_dbs{tri});
    end
    clear lfptrial_m1 lfptrial_dbs
    
    
    % save
    save(fullfile(savefolder, savefilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'T_idxevent_lfp', 'fs_ma', 'T_idxevent_ma', ...
        'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    
    
    clear outputfoldername dateofexp tdtbk
    clear pdcond savefilename onedaypath
    clear lfptrial_cortical lfptrial_dbs fs_lfp T_idxevent_lfp
    clear fs_ma T_idxevent_ma smoothWspeed_trial Wpos_smooth_trial Wrist_smooth_trial T_dbsChn
    clear lfpdata
end

end

function [lfptrial_notDBS, lfptrial_dbs, fs_lfp, T_idxevent_lfp, fs_ma, T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial, T_dbsChn] = extract_lfptrial(onedaypath, tdtblock)
% extractlfptrial extract trials for LFP data
%
%  [lfptrial_cortical, lfptrial_dbs, chantbl_cortical, chantbl_dbs] =
%  extractlfptrial(onedaypath, block) return extracted LFP trials of
%  cortical/subcortical data, dbs data, cortical/subcortical channel
%  information and dbs channel information tables, only trials with both
%  good reach and return are returned
%
%  Example usage:
%   onedaypath = '/home/lingling/root2/Animals/Jo/Recording/Processed/DataDatabase/Jo_020416'
%   tdtblock = 3
%   [lfptrial_cortical, lfptrial_dbs, fs, idxeventtbl, chantbl_dbs] = extractlfptrial_(onedaypath, tdtblock);
%
%  Inputs:
%   onedaypath : one date folder
%   tdtblock: tdt block number
%
%  Used files:
%       lfpfile_lfp  -   .\LFP\Block-3\Jo_CR1_DT1_020216_Block-3_LFPch*.nex
%       lfpfile_dbs   -   .\DBSLFP\Block-3\Jo_CR1_DT1_020216_Block-3_DBSLFP.nex.nex
%       mafile        -   .\Block-3\Jo_20160202_3_cleaned_MA_SingleTargetKluver_Analyze2.mat
%
%  Outputs:
%        lfptrial_notDBS: lfp trials of not DBS channels
%                      cell {1, n_trial}, inside each cell chns * ntemp
%
%        lfptrial_dbs: lfp trials of dbs channels
%                      cell {1, n_trial}, inside each cell chns * ntemp, 1-8: STN, 9-16:GP
%        fs_lfp: sample rate for lfp data
%
%
%        T_idxevent_lfp: a table describes the index for events of target onset,
%                      reach onset, touch screen, return and mouth in the lfp trial
%
%        fs_ma: sample rate for ma data
%
%
%        T_idxevent_ma: a table describes the index for events of target onset,
%                      reach onset, touch screen, return and mouth in the ma data trial
%
%
%        smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial: MA trial data
%
%        T_dbsChn:  a table describes each dbs channel information
%
%  More Description:
%       trial length = max(each trial length) + t_bef + t_aft
%       t_bef: the time before target on (default: t_bef = 1)
%       t_aft: the time after mouth (default: t_aft = 0.5)
%       one trailis from 't_target - t_bef'  to 't_mouth - t_after'

%% add NexMatablFiles path

[~, codefolder, ~, ~] = exp_subfolders();
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

%% MA data
% read the MA data
mafolder = fullfile(onedaypath, ['Block-' num2str(tdtblock)]);
mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));


if length(mafilestruct) ~= 1
    
    disp([mafolder 'has ' num2str(length(mafilestruct)) ' files, skip!'])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs_lfp = []; T_idxevent_lfp = [];
    fs_ma = []; T_idxevent_ma = [];
    smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
    T_dbsChn = [];
    
    return
end

% load SingleTargetKluverMAData from  *SingleTargetKluver_Analyze2.mat
load(fullfile(mafolder, mafilestruct.name), 'SingleTargetKluverMAData');

% ma sample rate
fs_ma = SingleTargetKluverMAData.SR;

% time indices for target onset, reach onset, touch screen, return and mouth
TargetTimeix = SingleTargetKluverMAData.TargetTime;
ReachonsetTimeix = SingleTargetKluverMAData.ReachTimeix;
TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
MouthTimeix = SingleTargetKluverMAData.MouthTimeix;
[m, n] = size(TargetTimeix);
if m == 1 || n == 1
    TargetTimeix = reshape(TargetTimeix, [m * n, 1]);
end

% timeix_ma matrix: ntrials * 5,   index for event from SingleTargetKluverMAData
timeix_ma = [TargetTimeix ReachonsetTimeix TouchTimeix ReturnTimeix MouthTimeix];
clear TargetTimeix ReachonsetTimeix TouchTimeix ReturnTimeix MouthTimeix

% the tag of good reach trials
tag_goodreach = SingleTargetKluverMAData.goodix_reach;
% the tag of good return trials
tag_goodreturn = SingleTargetKluverMAData.goodix_return;

% extract indices of good trials (both have good reach and return)
idx_goodtrials = (tag_goodreach == 1) & (tag_goodreturn == 1);

% if no good trials can be found
if isempty(idx_goodtrials)
    disp(['no good trials are found in ' mafilestruct.name])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs_lfp = []; T_idxevent_lfp = [];
    fs_ma = []; T_idxevent_ma = [];
    smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
    T_dbsChn = [];
    
    return;
end

% only remain the good trials
timeix_ma = timeix_ma(idx_goodtrials, :);
clear idx_goodtrials tag_goodreach tag_goodreturn


% varNames and varTypes for all return T_idxevent_ma, T_idxevent_lfp
varNames_table = {'TargetTimeix', 'ReachonsetTimeix', 'TouchTimeix', 'ReturnTimeix', 'MouthTimeix'};
varTypes_table = {'double','double','double', 'double', 'double'};

% total n_trial and n_event
[n_trial, n_events] = size(timeix_ma);

% t_bef: time before target on, t_aft: time after mouth
t_bef = 1;
t_aft = 0.5;


%% extract T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial

idx_mastrs = timeix_ma(:, 1) - round(t_bef * fs_ma); % start ma index for each trial n_trial * 1
idx_maends = timeix_ma(:, 5) + round(t_aft * fs_ma); % end ma index for each trial n_trial * 1


% smoothWspeed_trial,  Wpos_smooth_trial, Wrist_smooth_trial:  1* n_trial cell, each cell n_temporal * 1(3) 
for triali = 1:n_trial
    smoothWspeed_trial{triali} = SingleTargetKluverMAData.smoothWspeed(idx_mastrs(triali) : idx_maends(triali));
    Wpos_smooth_trial{triali} = SingleTargetKluverMAData.Wpos_smooth(idx_mastrs(triali) : idx_maends(triali));
    Wrist_smooth_trial{triali} = SingleTargetKluverMAData.Wrist_smooth(idx_mastrs(triali) : idx_maends(triali), :);
end
T_idxevent_ma = table('Size', size(timeix_ma), 'VariableTypes', varTypes_table, 'VariableNames',varNames_table);
T_idxevent_ma{:, :} = timeix_ma - repmat(idx_mastrs, [1, n_events]);


%% LFP data
% read each channel data in  LFP data to lfpdata (initialized when is the 1st channel)
folder_cortical = fullfile(onedaypath, 'LFP', ['Block-' num2str(tdtblock)]);

% files in folder_cortical are empty
if isempty(dir(folder_cortical))
    disp([folder_cortical ' has no files!'])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs_lfp = []; T_idxevent_lfp = [];
    fs_ma = []; T_idxevent_ma = [];
    smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
    T_dbsChn = [];
    
    return;
end

filenames = extractfield(dir(folder_cortical), 'name');
match = cellfun(@(x)~isempty(regexp(x, ['LFPch[0-9]*.nex'], 'match')), filenames, 'UniformOutput', false); % match channel file
match = cell2mat(match);
nexnames = filenames(match);

% parse and sort the number of channel
for filei = 1:length(nexnames)
    filename = nexnames{filei};
    tmp = char(regexp(filename, 'ch[0-9]*+.nex', 'match'));
    chns(filei) = str2num(tmp(length('ch') + 1:end - length('.nex')));
    
    if filei == 1
        % get the prefix of the each nex file name
        tmp = split(filename, [num2str(chns(filei)) '.nex']);
        file_prefix = tmp{1};
    end
    
end

chns = sort(chns);
for ich = 1:length(chns)
    filename = [file_prefix num2str(chns(ich)) '.nex'];
    [nexlfp_cortical] = readNexFile(fullfile(folder_cortical, filename));
    
    % extract the number of the structure containing LFPchn* data
    name_list = extractfield(cell2mat(nexlfp_cortical.contvars), 'name');
    i_lfp = find(contains(name_list, 'LFP')); %% i.e nexlfp_utah.contvars name == 'LFPch1', or 'MUAch1'
    
    
    if ich == 1% first channel
        fs_lfp = nexlfp_cortical.contvars{i_lfp}.ADFrequency;
        
        % timeix_lfp: index for event from nexlfp_cortical
        timeix_lfp = round(timeix_ma / fs_ma * fs_lfp);
        idx_lfpstrs = timeix_lfp(:, 1) - round(t_bef * fs_lfp);
        idx_lfpends = timeix_lfp(:, 5) + round(t_aft * fs_lfp);
        

        T_idxevent_lfp = table('Size', size(timeix_lfp), 'VariableTypes', varTypes_table, 'VariableNames',varNames_table);
        T_idxevent_lfp{:, :} = timeix_lfp - repmat(idx_lfpstrs, [1, n_events]);
        
        
        % lfptrial_cortical initial
        lfptrial_notDBS =  cell(1, n_trial);
    else     
        if fs_lfp ~= nexlfp_cortical.contvars{i_lfp}.ADFrequency% samping frequency is different
            chni = chns(ich);
            disp(['sampling frequency is different for chni = ' num2str(chni)]);
            break;
        end
    end
    
    
    
    % extract each trial for lfp data stored in separate channel
    for triali = 1:n_trial
        
        if size(nexlfp_cortical.contvars{i_lfp}.data, 1) < idx_lfpends(triali)
            disp(['chn = ' num2str(chns(ich)) ' LPF length < idx_end in triali = ' num2str(triali)])
            continue;
        end
        
        tmp = nexlfp_cortical.contvars{i_lfp}.data(idx_lfpstrs(triali):idx_lfpends(triali));
        tmp = reshape(tmp, 1, length(tmp));
        lfptrial_notDBS{triali} = [lfptrial_notDBS{triali}; tmp];   
        clear tmp
    end
    
    clear filename nexlfp_cortical name_list i_lfp
end



%% DBSLFP data
% read DBSLFP in nex file
dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtblock)]); % 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\DBSLFP\Block-1'

dbsfilepattern = fullfile(dbslfpfolder, '*DBSLFP.nex');
dbsfiles = dir(dbsfilepattern);

% dbs file does not exist
if length(dbsfiles) ~= 1
    disp([dbsfilepattern ' has ' num2str(length(dbsfiles)) ' file, skip!'])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs_lfp = []; T_idxevent_lfp = [];
    fs_ma = []; T_idxevent_ma = [];
    smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
    T_dbsChn = [];
    
    return;
end

[nexlfp_dbs] = readNexFile(fullfile(dbslfpfolder, dbsfiles(1).name));

% find the STN and GP data in the struct nexlfp_dbs
% extract all the values of field 'name' in nexlfp_dbs.convars
convars = cell2mat(nexlfp_dbs.contvars);
name_list = extractfield(cell2mat(nexlfp_dbs.contvars), 'name');
% extract the indices representing STN and GP data in nexlfp_dbs.contvars (33*1 struct array)
idx_dbs = find(contains(name_list, 'RAW_DBSch'));

% check the dbs channel frequencies
if range(cell2mat({convars(idx_dbs).ADFrequency})) ~= 0
    % check the sampling frequencies in channels are consistent
    disp(['ADFrequency in the dbs channels (STN & GP) is not consistent']);
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs_lfp = []; T_idxevent_lfp = [];
    fs_ma = []; T_idxevent_ma = [];
    smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
    T_dbsChn = [];
    
    return
end

if ~(abs(convars(idx_dbs(1)).ADFrequency - fs_lfp) < 0.001)
    % check the sampling frequency of dbs is the same as the lfp stored in separate channels or not
    disp(['the sampling frequency of dbs is not the same as the lfp stored in separate channels'])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs_lfp = []; T_idxevent_lfp = [];
    fs_ma = []; T_idxevent_ma = [];
    smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
    T_dbsChn = [];
    
    return
end


% extract each trial for lfp dbs data
nchn_dbs = length(idx_dbs);
lfptrial_dbs =  cell(1, n_trial);

for ich = 1:nchn_dbs
    chni = idx_dbs(ich);
    lfp_1chn = convars(chni).data;
    
    for triali = 1:n_trial
        tmp = lfp_1chn(idx_lfpstrs(triali):idx_lfpends(triali));
        tmp = reshape(tmp, 1, length(tmp));
        lfptrial_dbs{triali} = [lfptrial_dbs{triali}; tmp];  
        clear tmp
    end
    
    clear chni lfp_1chn triali
end

% chantbl_dbs: a table describes each dbs channel information
elecchn = extractfield(convars(idx_dbs), 'name');
elecchn = elecchn';
area = cell(nchn_dbs, 1);
area(1:8) = {'STN'};
area(9:16) = {'GP'};
T_dbsChn = [table(area) table(elecchn)];
clear varName

clear convars filename idx_stn idx_gp
end


function [T_chnsarea] = chanInf_M1DBS(nM1, nSTN, nGP)
% add M1 and DBS channel together
%
%   Args:
%       nM1: the number of M1 channels
%       nstn, ngp: the number of stn and gp channels
%
%   Output:
%       T_chnsarea: table of M1 + DBS channel inf,
%                   (T_chnsarea.Properties.VariableNames:
%                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

% initial the attitudes of chanarea table T_chnsarea:
% chni_vec, electype, brainareas, notes, recordingchn
chni_vec = uint8([1:nM1]');
electype = cell(nM1, 1);
brainareas = cell(nM1, 1);
notes = cell(nM1, 1);
recordingchn = uint8([1:nM1]');
electype(1:nM1, 1) = {'Gray Matter'}; % electype
brainareas(1:nM1, 1) = {'M1'}; % electype

% channel information table of M1
T_chnsarea = table;
T_chnsarea.chni = chni_vec;
T_chnsarea.brainarea = brainareas;
T_chnsarea.recordingchn = recordingchn;
T_chnsarea.electype = electype;
T_chnsarea.notes = notes;

% add DBS lead
for chi = 1:nSTN
    T_chnsarea = [T_chnsarea; {nM1 + chi, 'STN', chi, 'DBS', ''}];
end

for chi = 1:nGP
    T_chnsarea = [T_chnsarea; {nM1 + nSTN + chi, 'GP', chi, 'DBS', ''}];
end

end
