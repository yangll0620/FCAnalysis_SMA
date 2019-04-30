function extractlfptrial_m1dbs(animal, dateofexp, block, condition, therapy, badtrials)
% stop at the therapy part, Apr.26 2019, should identify this next week
%
%% extract trials for m1 and dbs lfp
%  now just can be used in Windows
%    
%  Example usage: 
%          extractlfptrial_m1dbs('Jo', '10-13-15', 8 ,'normal')      
%
%  Inputs:
%   dateofexp: format ('10-13-15', 'mm-dd-yy')
%   condition: 'normal', 'mild' or 'moderate'
%   therapy:   ''
%   badtrials: the trials will not be used.
%   load files:
%       lfpfile_utah  -   Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\LFP\Block-8\Jo_CR1_DT1_101315_Block-8_LFPch*.nex 
%       lfpfile_dbs   -   Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\DBSLFP\Block-8\Jo_CR1_DT1_101315_Block-8_DBSLFP.nex
%       mafile        -   Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\Block-8\Jo20151013_8_cleaned_MA_SingleTargetKluver_Analyze2.mat
%               
%
%  saveto: 
%       H:\My Drive\NMRC_umn\Projects\FunctionalConnectivityAnalysis\data\Jo\
%  savefile: 
%       Jo_lfptrial_m1dbs_normal_101315_block8.mat
%                                         chantbl_dbs (14x2 table):  describes area ('STN' or 'GP') and elecchn ('RAW_DBS_STNch3-4')
%             lfptrial_dbs (chn_dbs * n_temporal * n_trial double):  lfp trial data for dbs channels              
%           lfptrial_utah (chn_utah * n_temporal * n_trial double):  lfp trial data for utah data (for jo, all in M1) 
%                                     idxtbl_lfptrial (80x5 table):  stores the idx for events of target onset, reach onset, touch screen,
%                                                                    return and mouth in the trial matrix lfpdata_utah 
%                                                                    the first sample corresponds to target onset - t_bef               
%  More Description:
%       trial length = max(each trial length) + t_bef + t_aft
%       t_bef: the time before target on (default: t_bef = 1)
%       t_aft: the time after mouth (default: t_aft = 0.5)
%       one trailis from 't_target - t_bef'  to 't_mouth - t_after'


%% start from here
if strcmp(condition, 'normal') % no therapy for normal condition
    therapy = [];
end


%% file path setup
chn_lfputah = 96;
if ispc
    drivedir = 'Z:';
end
if isunix
    drivedir = '/run/user/1000/gvfs/ftp:host=nmrc_dserver1.local/root';
end
datapath = fullfile(drivedir, 'Animals', animal, 'Recording', 'Processed', 'DataDatabase'); % Z:\Animals\Jo\Recording\Processed\DataDatabase
onedaypath = fullfile(datapath, [animal '_' datestr(dateofexp,'mmddyy')]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_111715

%% add NexMatablFiles path
addpath(genpath(fullfile(pwd, 'toolbox', 'NexMatlabFiles')))


%% MA data
% read the MA data
mafolder = fullfile(onedaypath, ['Block-' num2str(block)]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\Block-8
mafile = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
%filename = [animal datestr(dateofexp, 'yyyymmdd') '_' num2str(block) '_cleaned_MA_SingleTargetKluver_Analyze2.mat'];
load(fullfile(mafolder, mafile.name)) % load SingleTargetKluverMAData
fs_ma = SingleTargetKluverMAData.SR;

% time indices for target onset, reach onset, touch screen, return and mouth
TargetTime = SingleTargetKluverMAData.TargetTime;
ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
MouthTimeix = SingleTargetKluverMAData.MouthTimeix;

timeixtbl_ma = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];
timeixtbl_ma(badtrials,:) = []; % remove the bad reaching trials
clear TargetTime ReachTimeix TouchTimeix ReturnTimeix MouthTimeix

%% utah array  LFP data
% folder_lfputah: utah lfp data path ( i.e. Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\LFP\Block-8)
folder_lfputah = fullfile(onedaypath, 'LFP', ['Block-' num2str(block)]);


t_bef = 1; t_aft = 0.5; % t_bef: time before target on, t_aft: time after mouth
[n_trial, n_timevars] = size(timeixtbl_ma);

% read each channel data in utah array  LFP data to lfpdata_utah (initialized when is the 1st channel)
for chni = 1 : chn_lfputah
    filename = [animal '_CR1_DT1_' datestr(dateofexp, 'mmddyy') '_Block-' num2str(block) '_LFPch' num2str(chni) '.nex']; % Jo_CR1_DT1_101315_Block-8_LFPch1.nex
    [nexlfp_utah] = readNexFile(fullfile(folder_lfputah, filename));
    
    % extract the number of the structure containing LFPchn* data
    name_list = extractfield(cell2mat(nexlfp_utah.contvars),'name');
    i_lfp = find(contains(name_list, 'LFP')); % % i.e nexlfp_utah.contvars name == 'LFPch1', or 'MUAch1'
    
    if chni == 1 % first channel
        fs_lfputah = nexlfp_utah.contvars{i_lfp}.ADFrequency;
        
        % the time index in the LFP neural data based on MA data
        timeixtbl_lfputah = timeixtbl_ma;
        timeixtbl_lfputah{:,:} = round(timeixtbl_ma{:,:} / fs_ma * fs_lfputah);
        
        % initialize lfp_utah : chn_lfputah * n_temporal * n_trial
        maxlen = max(timeixtbl_lfputah.MouthTimeix - timeixtbl_lfputah.TargetTime) + 1; % maximum length across all trials (unit: ind)
        n_bef = round(t_bef * fs_lfputah);  % n_bef: index number before target on
        n_aft = round(t_aft * fs_lfputah);  % n_aft: index number after mouth
        n_temporal = maxlen + n_bef + n_aft;
        idx_str  = timeixtbl_lfputah.TargetTime - n_bef;
        idx_end  = timeixtbl_lfputah.TargetTime + (maxlen -1) + n_aft;
        lfptrial_utah = zeros(chn_lfputah, n_temporal, n_trial);
        
        % idxtbl_lfptrialutah: the idx for events of target onset, reach onset, touch screen,
        % return and mouth in the trial matrix lfpdata_utah (chn_lfputah * n_temporal * n_trial)
        % the first sample corresponds to target onset - t_bef
        idxtbl_lfptrialutah = timeixtbl_lfputah;
        idxtbl_lfptrialutah{:,:} = idxtbl_lfptrialutah{:,:} - repmat(idx_str,[1, n_timevars]);
        
        clear timeixtbl_lfputah n_bef n_aft 
    else
        if fs_lfputah ~= nexlfp_utah.contvars{i_lfp}.ADFrequency % samping frequency is different
            disp(['sampling frequency is different for chni = ' num2str(chni)]);
            break;
        end
    end
    
    % extract each trial for lfp utah data
    for triali = 1: n_trial
        if size(nexlfp_utah.contvars{i_lfp}.data,1) < idx_end(triali)
            disp(mafile)
            disp(filename)
            disp(fs_ma)
            disp(fs_lfputah)
            disp(['idx in ma' num2str(timeixtbl_ma{1,end})])
            disp(size(nexlfp_utah.contvars{i_lfp}.data,1))
            disp(idx_end(triali))
            disp(['maxlen  = ' num2str(maxlen)])
        end
        lfptrial_utah(chni, :, triali) = nexlfp_utah.contvars{i_lfp}.data(idx_str(triali):idx_end(triali));
        
    end
    
    clear filename i_lfp 
end

% disp play the max trial number
disp(folder_lfputah)
disp(['max trial time is ' num2str(maxlen/fs_lfputah)]);

%% DBSLFP data
% read DBSLFP in nex file
dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(block)]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\DBSLFP\Block-8
filename = [animal '_CR1_DT1_' datestr(dateofexp, 'mmddyy') '_Block-' num2str(block) '_DBSLFP.nex']; % Jo_CR1_DT1_101315_Block-8_DBSLFP.nex
[nexlfp_dbs] = readNexFile(fullfile(dbslfpfolder, filename));

% find the STN and GP data in the struct nexlfp_dbs
% extract all the values of field 'name' in nexlfp_dbs.convars
convars = cell2mat(nexlfp_dbs.contvars);
name_list = extractfield(cell2mat(nexlfp_dbs.contvars),'name');
% extract the indices representing STN and GP data in nexlfp_dbs.contvars (33×1 struct array)
idx_stn = find(contains(name_list, 'RAW_DBS_STNch')); %
idx_gp = find(contains(name_list, 'RAW_DBS_GPch'));
idx_dbs = [idx_stn idx_gp];

% check the dbs channel frequencies
if range(cell2mat({convars(idx_dbs).ADFrequency})) ~= 0
    % check the sampling frequencies in channels are consistent
    disp(['ADFrequency in the dbs channels (STN & GP) is not consistent']);
    return
end
if(convars(idx_dbs(1)).ADFrequency ~= fs_lfputah) 
    % check the sampling frequency of dbs is the same as the lfp utah data or not
    disp(['the sampling frequency of dbs is not the same as the lfp utah data'])
    return
end
fs_lfpdbs = fs_lfputah; % sampling frequency for dbs channels
idxtbl_lfptrial = idxtbl_lfptrialutah;
clear idxtbl_lfputah

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

% chantbl_dbs: a table describes each channel information
elecchn = extractfield(convars(idx_dbs),'name');
elecchn = elecchn';
area = cell(nchn_dbs, 1);
area(1:length(idx_stn)) = {'STN'};
area(length(idx_stn)+1: end) = {'GP'};
chantbl_dbs =[table(area) table(elecchn)];
clear varName

clear convars filename idx_stn idx_gp 

if ispc % windows
end

savedir = 'H:\My Drive\NMRC_umn\Projects\FunctionalConnectivityAnalysis\data';
savepath = fullfile(savedir, animal);
savefile = [animal '_lfptrial_m1dbs_' condition therapy '_' datestr(dateofexp,'mmddyy') '_block' num2str(block) '.mat'];
save(fullfile(savepath, savefile), 'lfptrial_utah', 'lfptrial_dbs', 'idxtbl_lfptrial', 'chantbl_dbs')
