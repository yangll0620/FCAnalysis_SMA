function m0_SKTData_extract()
%% extract the STK data
%
%	Inputs:
%
%		xlsxfile_master, dateBlocks_w8-16UsedM1PMC.mat
%
%   Steps:
%       1. extract trials from processed folder in server (call function _extractlfptrial
%           a. trial length = max(each trial length) + t_bef + t_aft
%                   t_bef: the time before target on (default: t_bef = 1)
%                   t_aft: the time after mouth (default: t_aft = 0.5)
%           b. only remain the trials markedin both good reach and good return
%
%       2. Chns whose brain area are empty in Gray Matter are removed
%          For M1 and PMC, only chns within [8 16] are kept.
%
%       3. bipolar for DBS channels
%
%       4. Down sample trials into fs_new = 500
%
%       5. add variable T_chnsarea
%           




%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% codefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% datafolder, pipelinefolder 
[datafolder, ~, ~, ~] = exp_subfolders();
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

% channel number of grayMatter, DBS
nSTN = 8; nGP = 8; nGM = 96;

% Input dir:  preprocessed folder in root2
inputfolder_normalMild = fullfile(datafolder,animal, 'Recording', 'Processed', 'DataDatabase');
inputfolder_moderate = fullfile('Y:', 'Animals3', animal, 'Recording', 'Processed', 'DataDatabase');


% master sheet
xlsxfile_master = fullfile(datafolder, animal, [animal 'MasterDatabase.xlsx']);
% setup for normal
sheet_normal = 'Depth of GM array_Normal Channe';
nCols_normal = 36; % the total column number
row_area_normal = 19; row_chn_normal = 20;   % the row number for recording area, and channel number
% setup for mild
sheet_mild = 'Depth of GM array_Mild Channels';
nCols_mild = 93; % the total column number
row_area_mild = 2;     row_chn_mild = 3;   % the row number for recording area, and channel number
% setup for moderate
sheet_moderate = 'Depth of GM array_Moderate Chan';
nCols_moderate = 94; % the total column number
row_area_moderate = 1;          row_chn_moderate = 2;   % the row number for recording area, and channel number


% t_dateBlocks_SKT
file_dateBlocks = fullfile(codecorresParentfolder, '..', 'm2_dateBlocks_wUsedM1PMCChns', 'dateBlocks_w8-16UsedM1PMC.mat');
strformat_date_dateBlocks = 'yyyymmdd'; % the format of date string in master sheet, e.g '012317'


% downsample
fs_new = 500;

% conds
conds_cell = cond_cell_extract(animal);

%% save setup
savefolder = codecorresfolder;
savefilename_prefix = [animal '_STKData_'];


strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'


%% Starting

% ---  extract the same T_chnsarea_GM_noDepth and T_chnsarea_DBS for each day ----%

% T_chnsarea_GM
[T_chnsarea_GM_normal]  = chanInf_GM(xlsxfile_master, sheet_normal, nCols_normal, row_chn_normal, row_area_normal, nGM);
[T_chnsarea_GM_mild]  = chanInf_GM(xlsxfile_master, sheet_mild, nCols_mild, row_chn_mild, row_area_mild, nGM);
[T_chnsarea_GM_moderate]  = chanInf_GM(xlsxfile_master, sheet_moderate, nCols_moderate, row_chn_moderate, row_area_moderate, nGM);


% T_chnsarea_DBS for bipolar DBS
T_chnsarea_DBS= chanInf_DBS(nSTN - 1, nGP -1); 



% ---- extract t_SKT containing the records for STK ---- %
% extract t_dateBlocks_SKT_noDBS
load(file_dateBlocks, 't_dateBlocks_SKT_noDBS');



% ---- extract all STK trials ---- %
close all
f = waitbar(0, ['Extracting all STK trials']);
nrecords = height(t_dateBlocks_SKT_noDBS);
for i = 1 : nrecords
    % waitbar
    waitbar(i/nrecords,f,['i = ' num2str(i) ', Extracting trials in file ' num2str(i) '/' num2str(nrecords)]);
    
    
    %%% meta data for current record %%%%
    
    % date of exp, bktdt
    dateofexp = datenum(t_dateBlocks_SKT_noDBS.Date(i), strformat_date_dateBlocks);
    tdtbk = t_dateBlocks_SKT_noDBS.SKT_tdtbk(i);
    chnsUsed_M1 = t_dateBlocks_SKT_noDBS.chnsUsed_M1{i};
    chnsUsed_PMC = t_dateBlocks_SKT_noDBS.chnsUsed_PMC{i};
    
    
    % get the pd conditioon for the date of experiment
    pdcond = parsePDCondition(dateofexp, animal);
    if ~strcmp(pdcond, 'moderate')
        continue;
    end
    

    if ~any(strcmp(pdcond, conds_cell)) % avoid tomild, tomoderate
        continue;
    end
    
    
    
    %%% extract oneday+tdtbk data %%%
    savefilename = [savefilename_prefix  pdcond '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk) '.mat'];
    savefile = fullfile(savefolder, savefilename);

    
    switch pdcond
        case 'normal'
            inputfolder = inputfolder_normalMild;
        case 'mild'
            inputfolder = inputfolder_normalMild;
        case 'moderate'
            inputfolder = inputfolder_moderate;
        otherwise
            continue;
    end
        
    % one day path
    onedaypath = fullfile(inputfolder, [animal '_' datestr(dateofexp, 'mmddyy')]);
    
    if ~exist(onedaypath, 'dir')
        clear dateofexp tdtbk chnsUsed_M1 chnsUsed_PMC pdcondition
        clear onedaypath
        continue;
    end
    
    if exist(savefile, 'file') % if already dealed
        clear dateofexp tdtbk chnsUsed_M1 chnsUsed_PMC pdcond
        clear savefilename savefile
        continue;
    end
    
    % check if mafile *_Analyze2.mat/ LFP/DBSLFP all exist
    mafolder = fullfile(onedaypath, ['Block-' num2str(tdtbk)]); %  'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\Block-1'
    mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
    folder_cortical = fullfile(onedaypath, 'LFP', ['Block-' num2str(tdtbk)]);
    dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtbk)]);
    dbsfilepattern = fullfile(dbslfpfolder, ['*Block-' num2str(tdtbk) '_DBSLFP.nex']);
    if isempty(mafilestruct) || isempty(dir(folder_cortical))|| isempty(dir(dbsfilepattern))
        continue;
    end
    
    disp([onedaypath '-tdtbk' num2str(tdtbk)])
    continue;
    
    % extract trials of lfp data for particular date and tdt block number, lfptrial_*: nchns * ntemp * ntrials
    [lfptrial_GM, lfptrial_dbs,fs,T_idxevent, T_dbsChn] = extractlfptrial_(onedaypath, tdtbk);
    
    % skip this day if any is empty 
    if isempty(lfptrial_GM) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) ||isempty(T_dbsChn)
        clear dateofexp tdtbk chnsUsed_M1 chnsUsed_PMC pdcondition
        clear onedaypath
        clear lfptrial_GM lfptrial_dbs fs T_idxevent T_dbsChn
        
        continue;
    end
    

    
    %%% extract the lfptrial_GM used for brain areas %%%
    lfptrial_GM = lfptrial_GM([1:nGM], :,:);
    
    if(strcmp(pdcond, 'normal'))
        T_chnsarea_GM = T_chnsarea_GM_normal;
    else
        if (strcmp(pdcond, 'mild'))
            T_chnsarea_GM = T_chnsarea_GM_mild;
        else
            if (strcmp(pdcond, 'moderate'))
                T_chnsarea_GM = T_chnsarea_GM_moderate;
            end
        end
    end

    % the rows with brainarea is empty
    nonempty_mask = cellfun(@(x) ~isempty(x), T_chnsarea_GM.brainarea);
    % notusedM1_mask and  notusedPMC_mask
    notusedM1_mask = cellfun(@(x) strcmp(x, 'M1'), T_chnsarea_GM.brainarea) & ~ismember(T_chnsarea_GM.recordingchn, chnsUsed_M1);
    notusedPMC_mask = cellfun(@(x) strcmp(x, 'PMC'), T_chnsarea_GM.brainarea) & ~ismember(T_chnsarea_GM.recordingchn, chnsUsed_PMC);
 
    
    T_chnsarea_GM = T_chnsarea_GM(nonempty_mask & ~notusedM1_mask & ~notusedPMC_mask, :);
    lfptrial_GM = lfptrial_GM(nonempty_mask & ~notusedM1_mask & ~notusedPMC_mask, :, :);
    clear nonempty_mask notusedM1_mask notusedPMC_mask
    
    %%% bipolar dbs %%%%
    lfptrial_stn = diff(lfptrial_dbs(1:nSTN, :, :), [], 1);
    lfptrial_gp = diff(lfptrial_dbs(nSTN+1:end, :, :), [], 1);
    lfptrial_dbs = cat(1, lfptrial_stn, lfptrial_gp);
    clear lfptrial_stn lfptrial_gp
   
    
    
    %%% down sample %%%
    
    % lfptrial_GM down sample
    nchns = size(lfptrial_GM, 1);
    for chi = 1: nchns
        tmp = squeeze(lfptrial_GM(chi, :, :));
        tmp = resample(tmp, round(fs_new), round(fs));
        
        if chi == 1
            [ntemp, ntrials] = size(tmp);
            lfpdown_GM = zeros(nchns, ntemp, ntrials);
            
            clear ntemp ntrials
        end
        lfpdown_GM(chi, :, :) = tmp;
        
        clear tmp
    end
    
    % lfptrial_dbs down sample
    nchns = size(lfptrial_dbs, 1);
    for chi = 1: nchns
        tmp = squeeze(lfptrial_dbs(chi, :, :));
        tmp = resample(tmp, round(fs_new), round(fs));
        
        if chi == 1
            [ntemp, ntrials] = size(tmp);
            lfpdown_dbs = zeros(nchns, ntemp, ntrials);
            
            clear ntemp ntrials
        end
        lfpdown_dbs(chi, :, :) = tmp;
        
        clear tmp
    end
    
    % T_idxevent
    T_idxevent{:, :} =  round(T_idxevent{:, :} * fs_new / fs);
    fs = fs_new;
    

    %%% concatenate the Gray Matter chns, and the dbs chns  %%%
    T_chnsarea = vertcat(T_chnsarea_GM,T_chnsarea_DBS); T_chnsarea.chni = [1:height(T_chnsarea)]';
    lfpdata = cat(1, lfpdown_GM, lfpdown_dbs);
    clear T_chnsarea_GM lfpdown_dbs
    
    
    
    %%%  save  %%%
    save(savefile, 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    
    
    clear dateofexp tdtbk chnsUsed_M1 chnsUsed_PMC pdcondition
    clear onedaypath inputfolder
    clear lfptrial_GM lfptrial_dbs fs T_idxevent T_dbsChn
    clear lfpdata T_chnsarea savefilename savefile
   
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
    disp([mafolder ' has no Analyze2.mat'])
    
    lfptrial_cortical = []; lfptrial_dbs = [];
    fs = [];
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end


% load SingleTargetKluverMAData
load(fullfile(mafolder, mafilestruct.name), 'SingleTargetKluverMAData'); 

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
disp(['max trial time for today is ' num2str(maxlen/fs_lfpcortical)]);



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


function [T_chnsarea]  = chanInf_GM(xlsfile, sheetname, nCols, row_chn, row_area, nGM)
%
%   Extract T_chnsarea of Gray Matter for Bug
%
%   Inputs
%       xlsfile: the xls file (e.g /home/Bug/BugMasterDatabase.xlsx)
%       sheetname: the sheet name (e.g Depth of GM array_Moderate Chan)
%       nCols: the total column number, a scalar
%       row_chn:  the row number for chn
%       row_area: the row number for area
%       nGM: the total channel number of Gray Matter
%
%
%   Return:
%       t_areaDailyDepth: a table containing both area and daily depth (can be written directy using writetable)
%                         the variableNames are the date and channel numbers (e.g Date, chan4, chan19 et al)


% extract chnNames and areaNames
chnNames = readcell(xlsfile, 'FileType', 'spreadsheet', 'Sheet', sheetname, 'Range', ['A' num2str(row_chn) ':' colNum2ExcelColName(nCols) num2str(row_chn)]);
areaNames = readcell(xlsfile, 'FileType', 'spreadsheet', 'Sheet', sheetname, 'Range', ['A' num2str(row_area) ':' colNum2ExcelColName(nCols) num2str(row_area)]);
for i = 2: length(areaNames)
    if ismissing(areaNames{i})
        areaNames{i} = areaNames{i -1};
    end
end
areaNames{1} = '';

if(length(areaNames) ~= length(chnNames))
    disp('length of areaNames and chnNames not equal!')
    T_chnsarea =  [];
    return
end

% correspond chn with area
brainarea(1:nGM, 1) = {''};
for i = 2: length(chnNames)
    if isa(chnNames{i},'char')
        tmp = regexp(chnNames{i},'[0-9]*', 'match');
        if length(tmp) ~= 1
            disp([chnNames{i} ' matches not only 1 number'])
            T_chnsarea =  [];
            return
        end
        chni = str2num(tmp{1});
    end
    if isa(chnNames{i},'double')
        chni = chnNames{i};
    end
    
    if(~isempty(brainarea{chni,1}))
        disp([num2str(chni) ' is extracted twice: ' areaNames{i} ' and ' brainarea{chni,1}])
        T_chnsarea =  [];
        return
    end
    
    brainarea(chni,1) = areaNames(i);
end 


% generate T_chnsarea
T_chnsarea = table;
T_chnsarea.chni = uint8([1:nGM]');
T_chnsarea.brainarea = brainarea;
T_chnsarea.recordingchn = T_chnsarea.chni;
T_chnsarea.electype(1:nGM, 1) = {'Gray Matter'};
T_chnsarea.notes(1:nGM, 1) = {''};

end % end chanInf_GM


function cCol= colNum2ExcelColName(colNum)
% colnum to excel column name
% Exp: 1-> 'A',  27 -> 'AA'
% Input
%       colNum: the column Number
% output:
%       cCol: the transformed excel column name

div = floor(colNum/26);
re = mod(colNum, 26);


cRe = char(re + 'A' - 1);

if div ==  0
    cDiv = '';
else
    cDiv = char(div + 'A' - 1) ;
end
cCol = [cDiv cRe];

end %colNum2ExcelColName




function T_chnsarea = chanInf_DBS(nSTN, nGP)
    % extract M1 channel inf table
    %   Args:
    %       nSTN, nGP: the number of stn and gp channels
    %
    %   Output:
    %       T_chnsarea: table of DBS channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    % channel information table of M1
    T_chnsarea = table;

    T_chnsarea.chni = uint8([1:nSTN + nGP]');

    T_chnsarea.brainarea = cell(nSTN + nGP, 1);
    T_chnsarea.brainarea(1:nSTN) = {'STN'}; T_chnsarea.brainarea(nSTN + 1:nSTN + nGP) = {'GP'};

    T_chnsarea.recordingchn = [uint8([1:nSTN]'); uint8([1:nGP]')];

    T_chnsarea.electype = cell(nSTN + nGP, 1);
    T_chnsarea.electype(:) = {'DBS'};

    T_chnsarea.notes = cell(nSTN + nGP, 1);

end




