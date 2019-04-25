clear

%% input: animal, dateofexp, block et.al
chn = 96;
animal = 'Jo';
dateofexp = datenum('10-13-15', 'mm-dd-yy');
block = 8;
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

%% read utah array  LFP data locating on M1 in nex file
lfputahfolder = fullfile(onedaypath, 'LFP', ['Block-' num2str(block)]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\LFP\Block-8
for chni = 1 : 1%chn
    filename = [animal '_CR1_DT1_' datestr(dateofexp, 'mmddyy') '_Block-' num2str(block) '_LFPch' num2str(chni) '.nex']; % Jo_CR1_DT1_101315_Block-8_LFPch1.nex
    [nexlfp_utah] = readNexFile(fullfile(lfputahfolder, filename));
    
    for i = 1: length(nexlfp_utah.contvars)
        convname = nexlfp_utah.contvars{i}.name;
        if ~isempty(findstr(convname, 'LFP')) % i.e convname == 'LFPchn1', or convname == 'MUAch1'
            break;
        end
    end
    
    fs_lfputah = nexlfp_utah.contvars{i}.ADFrequency;
    lfp_utah = nexlfp_utah.contvars{i}.data; % .data: n_temporal * 1
    clear filename
end

%% read DBSLFP in nex file
dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(block)]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\DBSLFP\Block-8
filename = [animal '_CR1_DT1_' datestr(dateofexp, 'mmddyy') '_Block-' num2str(block) '_DBSLFP.nex']; % Jo_CR1_DT1_101315_Block-8_DBSLFP.nex
[nexlfp_dbs] = readNexFile(fullfile(dbslfpfolder, filename));
clear filename

%% read MA data
mafolder = fullfile(onedaypath, ['Block-' num2str(block)]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Jo_101315\Block-8
filename = [animal datestr(dateofexp, 'yyyymmdd') '_' num2str(block) '_cleaned_MA_SingleTargetKluver_Analyze2.mat'];
load(fullfile(mafolder, filename)) % load SingleTargetKluverMAData
fs_ma = SingleTargetKluverMAData.SR;

% time indices for target onset, reach onset, touch screen, return and mouth 
TargetTime = SingleTargetKluverMAData.TargetTime;
ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix 
MouthTimeix = SingleTargetKluverMAData.MouthTimeix;
matimeix_tbl = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];

%% calculate the time index in the neural data
neurtimeix_tbl = matimeix_tbl / fs_ma * fs_lfp;
neurtimeix_tbl = neurtimeix_tbl{:,2:end} = T{:,2:end}*25/100



%% Extracting field values from cell array of structures - will use tmr, 4/25/2019
RP_Array = cell2mat(nexlfp_dbs.contvars);
fnames = fieldnames(RP_Array);
Allfields = struct2cell(RP_Array); 
name_List = squeeze(Allfields(strcmp('name', fieldnames(RP_Array)), :, :));


% code needed to be tested tmr, 4/25/2019
%  https://www.mathworks.com/matlabcentral/answers/295641-extracting-field-values-from-cell-array-of-structures
Start_List = cellfun(@(c) {c.Start}, RP_Cell, 'UniformOutput', false);

