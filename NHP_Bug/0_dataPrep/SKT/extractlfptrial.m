function [lfptrial_cortical, lfptrial_dbs, fs ,idxeventtbl, chantbl_dbs] = extractlfptrial(onedaypath, tdtblock)
% extractlfptrial extract trials for LFP data
%
%  [lfptrial_cortical, lfptrial_dbs, chantbl_cortical, chantbl_dbs] =
%  extractlfptrial(onedaypath, block) return extracted LFP trials of 
%  cortical/subcortical data, dbs data, cortical/subcortical channel
%  information and dbs channel information tables, only trials with both
%  good reach and return are returned
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
%       one trailis from 't_target - t_bef'  to 't_mouth - t_after'

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