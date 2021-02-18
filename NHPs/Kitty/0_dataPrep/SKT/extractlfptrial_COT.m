function [lfptrial_notDBS, lfptrial_dbs, fs, idxeventtbl, chantbl_dbs] = extractlfptrial_COT(onedaypath, tdtblock, file_SPK)
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
%       lfpfile_utah  -   .\LFP\Block-3\Jo_CR1_DT1_020216_Block-3_LFPch*.nex
%       lfpfile_dbs   -   .\DBSLFP\Block-3\Jo_CR1_DT1_020216_Block-3_DBSLFP.nex.nex
%       mafile        -   .\Block-3\Jo_20160202_3_cleaned_MA_SingleTargetKluver_Analyze2.mat
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




%%%--- extract MA information from SPK nex file ---%%%
events_COT = {'TargetAppear', 'ReachStart', 'Touch', 'ReturnStart', 'Mouth'};
id_target = 1;

nex_SPK = readNexFile(file_SPK);

tbeg_SPK = nex_SPK.tbeg;
tend_SPK = nex_SPK.tend;
events_SPK = nex_SPK.events;

for ei = 1: length(events_SPK)
    event = events_SPK{ei}.name;
    event = strrep(event, num2str(id_target), '');
    
    if any(strcmp(event, events_COT))
        eval(['time_' event ' = events_SPK{ei}.timestamps;'])
    end
    
    clear event
end


timetbl_SPK = table;
for evi = 1 : length(events_COT)
    event = events_COT{evi};
    eval(['timetbl_SPK = [timetbl_SPK table(time_' event ')];'])
    clear event
end





[n_trial, n_timevars] = size(timetbl_SPK);


%%%--- not DBS LFP data ---%%%
t_bef = 1; t_aft = 0.5; % t_bef: time before target on, t_aft: time after mouth

% read each channel data in  LFP data to lfpdata (initialized when is the 1st channel)
folder_notDBSLFP = fullfile(onedaypath, 'LFP', ['Block-' num2str(tdtblock)]);

% files in folder_cortical are empty
if isempty(dir(folder_notDBSLFP))
    disp([folder_notDBSLFP ' has no files!'])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs = [];
    idxeventtbl = []; chantbl_dbs = [];
    
    return;
end

filenames = extractfield(dir(folder_notDBSLFP), 'name');
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
chn_lfp = length(chns);


for chi = 1:length(chns)
    filename = [file_prefix num2str(chns(chi)) '.nex'];
    [nexlfp_notDBS] = readNexFile(fullfile(folder_notDBSLFP, filename));
    
    % extract the number of the structure containing LFPchn* data
    name_list = extractfield(cell2mat(nexlfp_notDBS.contvars), 'name');
    i_lfp = find(contains(name_list, 'LFP')); %% i.e nexlfp_utah.contvars name == 'LFPch1', or 'MUAch1'
    
    if chi == 1% first channel
        if tbeg_SPK ~= nexlfp_notDBS.tbeg || tend_SPK ~= nexlfp_notDBS.tend
            
            disp('tbeg or tend not consistent in SPK and LFP Nex files')
            
            lfptrial_notDBS = []; lfptrial_dbs = [];
            fs = [];
            idxeventtbl = []; chantbl_dbs = [];
            
            return;
        end
        
        
        fs_lfpNotDBS = nexlfp_notDBS.contvars{i_lfp}.ADFrequency;
        
        % the time index in the LFP neural data based on MA data
        timeixtbl_lfpNotDBS = timetbl_SPK;
        timeixtbl_lfpNotDBS{:, :} = round(timetbl_SPK{:, :} * fs_lfpNotDBS);
        timeixtbl_lfpNotDBS.Properties.VariableNames{'time_TargetAppear'} = 'TargetTimeix';
        timeixtbl_lfpNotDBS.Properties.VariableNames{'time_ReachStart'} = 'ReachTimeix';
        timeixtbl_lfpNotDBS.Properties.VariableNames{'time_Touch'} = 'TouchTimeix';
        timeixtbl_lfpNotDBS.Properties.VariableNames{'time_ReturnStart'} = 'ReturnTimeix';
        timeixtbl_lfpNotDBS.Properties.VariableNames{'time_Mouth'} = 'MouthTimeix';
        
        % initialize lfptrial_notDBS : chn_lfp * n_temporal * n_trial
        maxlen = max(timeixtbl_lfpNotDBS.MouthTimeix - timeixtbl_lfpNotDBS.TargetTimeix) + 1; % maximum length across all trials (unit: ind)
        n_bef = round(t_bef * fs_lfpNotDBS); % n_bef: index number before target on
        n_aft = round(t_aft * fs_lfpNotDBS); % n_aft: index number after mouth
        n_temporal = maxlen + n_bef + n_aft;
        idx_str = timeixtbl_lfpNotDBS.TargetTimeix - n_bef;
        idx_end = timeixtbl_lfpNotDBS.TargetTimeix + (maxlen -1) + n_aft;
        lfptrial_notDBS = zeros(chn_lfp, n_temporal, n_trial);
        
        % idxtbl_lfptrialutah: the idx for events of target onset, reach onset, touch screen,
        % return and mouth in the trial matrix lfpdata_utah (chn_lfputah * n_temporal * n_trial)
        % the first sample corresponds to target onset - t_bef
        idxtbl_lfptrial_notDBS = timeixtbl_lfpNotDBS;
        idxtbl_lfptrial_notDBS{:, :} = idxtbl_lfptrial_notDBS{:, :} - repmat(idx_str, [1, n_timevars]);
        
        clear timeixtbl_lfpsepchn n_bef n_aft
    else
        
        if fs_lfpNotDBS ~= nexlfp_notDBS.contvars{i_lfp}.ADFrequency% samping frequency is different
            chni = chns(chi);
            disp(['sampling frequency is different for chni = ' num2str(chni)]);
            break;
        end
        
    end
    
    % extract each trial for lfp data stored in separate channel
    for triali = 1:n_trial
        
        if size(nexlfp_notDBS.contvars{i_lfp}.data, 1) < idx_end(triali)
            disp(mafilestruct)
            disp(filename)
            disp(fs_ma)
            disp(fs_lfpNotDBS)
            disp(['idx in ma' num2str(timeixtbl_ma{1, end})])
            disp(size(nexlfp_notDBS.contvars{i_lfp}.data, 1))
            disp(idx_end(triali))
            disp(['maxlen  = ' num2str(maxlen)])
        end
        
        lfptrial_notDBS(chi, :, triali) = nexlfp_notDBS.contvars{i_lfp}.data(idx_str(triali):idx_end(triali));
        
    end
    
    clear filename i_lfp
end

%%%--- DBS LFP data ---%%%
% read DBSLFP in nex file
dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtblock)]); % 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\DBSLFP\Block-1'

dbsfilepattern = fullfile(dbslfpfolder, '*DBSLFP.nex');
dbsfiles = dir(dbsfilepattern);

% dbs file does not exist
if length(dbsfiles) ~= 1
    disp([dbsfilepattern ' has ' num2str(length(dbsfiles)) ' file, skip!'])
    
    lfptrial_notDBS = []; lfptrial_dbs = [];
    fs = [];
    idxeventtbl = []; chantbl_dbs = [];
    
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
    return
end

if (convars(idx_dbs(1)).ADFrequency ~= fs_lfpNotDBS)
    % check the sampling frequency of dbs is the same as the lfp stored in separate channels or not
    disp(['the sampling frequency of dbs is not the same as the lfp stored in separate channels'])
    return
end

% extract each trial for lfp dbs data
nchn_dbs = length(idx_dbs);
lfptrial_dbs = zeros(nchn_dbs, n_temporal, n_trial);

for chi = 1:nchn_dbs
    chni = idx_dbs(chi);
    lfp_1chn = convars(chni).data;
    
    for triali = 1:n_trial
        lfptrial_dbs(chi, :, triali) = lfp_1chn(idx_str(triali):idx_end(triali));
    end
    
    clear chni lfp_1chn triali
end

% chantbl_dbs: a table describes each dbs channel information
elecchn = extractfield(convars(idx_dbs), 'name');
elecchn = elecchn';
area = cell(nchn_dbs, 1);
area(1:8) = {'STN'};
area(9:16) = {'GP'};
chantbl_dbs = [table(area) table(elecchn)];
clear varName

clear convars filename idx_stn idx_gp


%%%---uniform fs and idxeventbl ---%%%
fs = fs_lfpNotDBS; % sampling frequency for dbs channels
idxeventtbl = idxtbl_lfptrial_notDBS;
clear idxtbl_lfptrial_notDBS


end

