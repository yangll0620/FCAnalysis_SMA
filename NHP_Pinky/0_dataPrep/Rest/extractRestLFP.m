function [lfpRest_UtahGM, lfpRest_DBS, fs, tbl_ChanDBS] = extractRestLFP(onedaypath, block)
% extract LFP rest data for all the UtahGrayMatter and DBS channels in oneday\block 
%
% Return:
%   lfpRest_UtahGM: nchns * ntemp
%
%   lfpRest_DBS: nchns * ntemp
%  
%   fs: sample rate
%
%   tbl_ChanDBS: dbs channel table

%% add NexMatablFiles path
addpath(genpath(fullfile('..', '..', 'toolbox', 'NexMatlabFiles')))


%% LFP data for Utah + GM channels
% nex folder for Utah+GM data 
nexfolder_UtahGM = fullfile(onedaypath, 'LFP', ['Block-' num2str(block)]);

[file_prefix, ~] = nexfilenames_infolder(nexfolder_UtahGM);
chns_UtahGM = [1:96 101:132];
for i = 1: length(chns_UtahGM)
    
    filename = [file_prefix num2str(chns_UtahGM(i)) '.nex'];
    
    % read nex file data
    [nexdata] = readNexFile(fullfile(nexfolder_UtahGM, filename));
    
    % extract the number of the structure containing LFPchn* data
    name_list = extractfield(cell2mat(nexdata.contvars),'name');
    
    % nexlfp.contvars name == 'LFPch1', or 'MUAch1'
    i_lfp = find(contains(name_list, 'LFP'));
    
    % lfp data
    data = nexdata.contvars{i_lfp}.data;
    
    if i == 1 
        % sample rate
        fs = nexdata.contvars{i_lfp}.ADFrequency;
        
        % lfpdata is for the lfp data of all the chns
        lfpRest_UtahGM = zeros(length(chns_UtahGM), length(data));
   
    else
        % samping frequency is different
        if fs ~= nexdata.contvars{i_lfp}.ADFrequency
            chni = chns_UtahGM(i);
            disp([nexfolder 'fs = ' num2str(fs) ',sampling frequency is different for chni = ' num2str(chni)]);
            lfpRest_UtahGM = [];
            fs = [];
            break;
        end
    end
    
    % assign to lfpdata
    lfpRest_UtahGM(i, :) = data;
    
    clear filename nexdata name_list i_lfp data
end
fs_UtahGM = fs;


%% DBS LFP data
% nex folder for DBS data 
nexfolder_DBS = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(block)]);

% dbs file does not exist
if isempty(dir(fullfile(nexfolder_DBS, '*_GrayMatter*DBSLFP.nex')))
    disp([dbslfpfolder, '*_GrayMatter*DBSLFP.nex files do not exit!'])
    
    lfpRest_UtahGM = []; lfpdata_DBS = [];
    fs = [];
    idxeventtbl = []; tbl_ChanDBS = [];
    
    return;
end

filenames = extractfield(dir(fullfile(nexfolder_DBS, '*_GrayMatter*DBSLFP.nex')),'name'); % Pinky_GrayMatter_eyetracking_DT1_071417_Block-1_DBSLFP.nex
match = cellfun(@(x) ~isempty(regexp(x, ['[A-Z, a-z, 0-9]*.nex'],'match')),filenames,'UniformOutput',false); % match channel file
match = cell2mat(match);
nexnames = filenames(match);
if length(nexnames) ~= 1
    % only has one file inside folder
    error(['More than one _DBSLFP.nex  in file' nexfolder_DBS])
end

% read DBS lfp data in nex file
[nexlfp_dbs] = readNexFile(fullfile(nexfolder_DBS, nexnames{1}));

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
if(convars(idx_dbs(1)).ADFrequency ~= fs_UtahGM) 
    % check the sampling frequency of dbs is the same as the lfp stored in separate channels or not
    disp(['the sampling frequency of dbs is not the same as the lfp stored in UtahGM channels'])
    return
end
fs = fs_UtahGM;

% extract each trial for lfp dbs data
nchn_dbs = length(idx_dbs);
for i = 1: nchn_dbs
    chni = idx_dbs(i);
    lfp_1chn = convars(chni).data;
    
    if i == 1
        n_temps = length(lfp_1chn);
        lfpRest_DBS = zeros(nchn_dbs, n_temps);
        clear n_temps
    end
    
    lfpRest_DBS(i, :) = lfp_1chn;    
    
    clear lfp_1chn chni 
end

% chantbl_dbs: a table describes each dbs channel information
elecchn = extractfield(convars(idx_dbs),'name');
elecchn = elecchn';
area = cell(nchn_dbs, 1);
area(1:8) = {'STN'};
area(9: 16) = {'GP'};
tbl_ChanDBS =[table(area) table(elecchn)];
clear varName

