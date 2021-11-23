function m1_restData_cleaned_extract_moderate()
%% extract cleaned restData using files Ying used
%
%   1. remove data of channels in m1array_to_remove
%   
%   2. remove data with eye-closed
%
%   3. remove data marked with movement 
%
%   4. Chns whose brain area are empty in Gray Matter are removed
%          For M1 and PMC, only chns within [8 16] are kept.
%   
%   Inputs:
%       
%       configFile: datafolder/ config_m1lf_fromYing.mat
%
%       all_sync file: e.g /Block-1/Bug_20180411_1_all_sync.mat



%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');


% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

% datafolder
[datafolder, ~, pipelinefolder, ~] = exp_subfolders();

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

% Input dir:  preprocessed folder in root2
inputfolder_normalMild = fullfile(datafolder,animal, 'Recording', 'Processed', 'DataDatabase');
inputfolder_moderate = fullfile('Y:', 'Animals3', animal, 'Recording', 'Processed', 'DataDatabase');
dateBlocks_file = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'm2_dateBlocks_wUsedM1PMCChns','dateBlocks_w8-16UsedM1PMC.mat');
% threshold used by Ying to extract cleaned resting data
configFile = fullfile(datafolder, 'root2', 'Ying Yu', 'SKB_EventRelatedBeta', 'PAC', 'config_m1lf.mat');


% channel number of grayMatter, DBS
nSTN = 8; nGP = 8; nGM = 96;

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


%% save setup
savefolder = codecorresfolder;
savefilename_prefix = [animal '_cleanedRestData_'];

%% global variables

% to be removed channels in array
arraychns_to_remove = [];


% start time point, data before t0 are removed
t0 = 3;

% only remain the interval whose duration is larger than min_comp_time
min_comp_time = 15;

% threshold for eye_open
thr_eye_open = 0.5000;


%%  Starting here

% ---  extract the same T_chnsarea_GM_noDepth and T_chnsarea_DBS for each day ----%

% T_chnsarea_GM
[T_chnsarea_GM_normal]  = chanInf_GM(xlsxfile_master, sheet_normal, nCols_normal, row_chn_normal, row_area_normal, nGM);
[T_chnsarea_GM_mild]  = chanInf_GM(xlsxfile_master, sheet_mild, nCols_mild, row_chn_mild, row_area_mild, nGM);
[T_chnsarea_GM_moderate]  = chanInf_GM(xlsxfile_master, sheet_moderate, nCols_moderate, row_chn_moderate, row_area_moderate, nGM);


% T_chnsarea_DBS for bipolar DBS
T_chnsarea_DBS= chanInf_DBS(nSTN, nGP); 



% ---- extract t_dateBlocks_Rest_noDBS containing the records for Rest ---- %
% load dateBlocks
load(dateBlocks_file, 't_dateBlocks_Rest_noDBS');

% load configFile
load(configFile, 'config_m1lf');


% ---- extract all Rest segments ---- %
nfiles = height(t_dateBlocks_Rest_noDBS);
for fi = 1: nfiles
    
    % Date and tdt block 
    dateofexp = datenum(t_dateBlocks_Rest_noDBS.Date(fi), 'yyyymmdd');
    tdtblock = t_dateBlocks_Rest_noDBS.Rest_tdtbk(fi);
    mablock = t_dateBlocks_Rest_noDBS.Rest_mabk(fi);
    pdcond = parsePDCondition(dateofexp, animal);
    
    switch pdcond
        case 'normal'
            continue;
        case 'mild'
            continue;
        case 'moderate'
            inputfolder = inputfolder_moderate;
        otherwise
            continue;
    end
    
    % sync file pattern
    folder_date = fullfile(inputfolder, [animal '_' datestr(dateofexp, 'mmddyy')]);
    if ~exist(folder_date, 'dir')
        continue;
    end
    
    
    folder_allsync = fullfile(folder_date,  ['Block-' num2str(tdtblock)]);
    filepattern_allsync = fullfile(folder_allsync, [animal '_' datestr(dateofexp, 'yyyymmdd') '_' num2str(mablock) '_all_sync.mat']);
    files_allsync = dir(filepattern_allsync);
    
    % all sync file not exist or there are more than one sync file, skip
    if(length(files_allsync) ~= 1)
        disp([filepattern_allsync ' has ' num2str(length(files_allsync)) ' files, not 1 file, skipped.'])
        continue;
    end
    
    % the all sync file 
    file_allsync = fullfile(files_allsync.folder, files_allsync.name);
    
    
    %%%  extract the segments %%%
        
    if strcmp(pdcond, 'mild') || strcmp(pdcond, 'moderate')
        eval(['thr_power = config_m1lf.max_' lower(animal) '_pd;'])
    else
        if strcmp(pdcond, 'normal')
            eval(['thr_power = config_m1lf.max_' lower(animal) '_normal;'])
        else
            disp([datestr(dateofexp, 'yyyymmdd') ', Ignore condition = ' pdcond])
            
            continue;
        end
    end
    
    disp([num2str(fi) '/' num2str(nfiles) ': extracting segments from ' datestr(dateofexp, 'yyyymmdd') '-Block' num2str(tdtblock)])
    [data_segments, segsIndex, therapy, fs] = segments_extract_fromAllSyncMat(file_allsync, arraychns_to_remove, thr_power, thr_eye_open, t0, min_comp_time, config_m1lf.freq_band);
    
    
    if(isempty(data_segments))
        continue;
    end
    
    
    %%% extract the data_segments used for brain areas %%%
    
    chnsUsed_M1 = t_dateBlocks_Rest_noDBS.chnsUsed_M1{fi};
    chnsUsed_PMC = t_dateBlocks_Rest_noDBS.chnsUsed_PMC{fi};
    
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
    chnsUsed_mask = nonempty_mask & ~notusedM1_mask & ~notusedPMC_mask;
    
    T_chnsarea_GM = T_chnsarea_GM(chnsUsed_mask, :); T_chnsarea_GM.chni = [1:height(T_chnsarea_GM)]';
    for segi = 1 : length(data_segments)
        data_segments(segi).lfp_array = data_segments(segi).lfp_array(:, chnsUsed_mask);
    end
    clear nonempty_mask notusedM1_mask notusedPMC_mask chnsUsed_mask chnsUsed_M1 chnsUsed_PMC
    
    
    
    
    str_therapy = '';
    if (strcmp(pdcond, 'mild') || strcmp(pdcond, 'moderate')) && ~strcmp(therapy, 'off')
            str_therapy = therapy;
    end
    
    savefile = fullfile(savefolder, [savefilename_prefix  pdcond str_therapy '_'  datestr(dateofexp,'yyyymmdd') '_tdt'  num2str(tdtblock) '.mat']);
    
    save(savefile, 'data_segments', 'segsIndex', 'fs', 'T_chnsarea_GM', 'T_chnsarea_DBS')
    
    clear allsyncfile dateblockstr thr_power1
    clear tmp condition 
    clear data_segments segsIndex fs T_chnsarea_GM 
end
end




function [data_segments, segsIndex, therapy, fs] = segments_extract_fromAllSyncMat(file_allsync, arraychns_to_remove, thr_power, thr_eye_open, t0, min_comp_time, freq_band)
% segments_extrct_fromAllSyncMat
%   extract the segments whose 
%
%
%   Args:
%       file_allsync: the all sync mat file (fullpath, e.g. /.../../Bug_20190110_2_all_sync.mat)
%       arraychns_to_remove: the array channels to be removed (e.g. [] or [1, 4])
%       thr_power1
%       thr_eye_open: 
%       t0
%       min_comp_time
%       freq_band
%
%
%   Returns:
%       data_segments:
%
%       segsIndex:
%   




%% Start here

% load all sync file
load(file_allsync, 'data');


fs = data.lfp_stn_fs(1);

if(fs < 3000)
    % fs < 3000, the designed filter in filtered_lfp() is not right
    
    data_segments = [];
    segsIndex = [];
    therapy = [];
end



n0 = round(t0 * fs);
min_samples = round(min_comp_time * fs);

% removing the arraychns_to_remove channels
lfp_array = zeros(length(data.GMdata{1}), length(data.GMdata));
for chi = 1 : length(data.GMdata)
    lfp_array(:, chi) = data.GMdata{chi};
end
lfp_array(:, arraychns_to_remove) = [];


% get lfp_stn, and lfp_gp data
lfp_stn = zeros(length(data.lfp_stn{1}),length(data.lfp_stn));
for chi =1:length(data.lfp_stn)
    lfp_stn(:,chi) = data.lfp_stn{chi};
end

lfp_gp = zeros( length(data.lfp_gp{1}),length(data.lfp_gp));
for chi =1:length(data.lfp_gp)
    lfp_gp(:,chi) = data.lfp_gp{chi};
end


%%  get states for each time point of lfp, 1: remain this time point, 0: remove
chns_m1 = [56 58 79 69 65];
avglfp_m1 = mean(lfp_array(:, chns_m1), 2);
m1power_state_mintime = extract_power_state_mintime(avglfp_m1, thr_power, min_samples, fs,freq_band);

ntemp = length(lfp_array);
time_us = [0:ntemp-1]/fs;
eye_state_mintime = extract_eye_state_mintime(data.eye_area_filt, thr_eye_open, time_us, min_samples);

% awake state should be both eye is open and m1_power is less than thr_power1
awake_state_mintime = ( (eye_state_mintime>=1) & (m1power_state_mintime>=1)) ;


% start and end index pair for each no movement segments idx_NoMoveSegs: n_NoMoveSegs * 2 (idx_start, idx_end)
idx_noMoveSegs = round(data.SponMovement.not_movement_times * fs);
to_remove_ma = extract_to_remove_ma(idx_noMoveSegs, ntemp);

states = awake_state_mintime & (~to_remove_ma);


%% filter lfp_m1, lfp_stn, lfp_gp, and truncate at n0
[lfp_array, lfp_stn, lfp_gp]= filtered_lfp(lfp_array, lfp_stn, lfp_gp, fs);


%% truncate
lfp_array(1: n0-1, :) = [];
lfp_stn(1: n0-1, :) = [];
lfp_gp(1: n0-1,:) = [];
states(1:n0-1) = [];


%% extract segments
%  resting data segments indices vec_segsIndex: nsegs * 2
segsIndex = get_segIndex( states, min_samples );

% extract segments
for seg_ind = 1: size(segsIndex,1)
    
    idx_segStr = segsIndex(seg_ind,1);
    idx_segEnd = segsIndex(seg_ind,2);
    
    %  lfp_array
    if ~exist('data_segments')
        data_segments(1).lfp_array = lfp_array(idx_segStr:idx_segEnd, :);
    else
        data_segments(end+1).lfp_array = lfp_array( idx_segStr:idx_segEnd, :);
    end
    
    % lfp_stn and lfp_gp
    data_segments(end).lfp_stn = lfp_stn(idx_segStr:idx_segEnd, :);
    data_segments(end).lfp_gp = lfp_gp(idx_segStr:idx_segEnd, :);
    
    
    clear idx_segStr idx_segEnd
end

if ~exist('data_segments')
    data_segments = [];
    segsIndex = [];
end

therapy = data.therapy;
end



function [lfp_array, lfp_stn, lfp_gp]= filtered_lfp(lfp_array, lfp_stn, lfp_gp, fs)
%% pass filter lfp_array, lfp_stn, lfp_gp  

% pass filters of array, dbs data
[bhp,ahp] = butter(2 , 2*2/fs , 'high'); % high pass


% high pass filter lfp_array
for chi = 1: size(lfp_array, 2)
    lfp_array(:,chi) = filtfilt(bhp,ahp, lfp_array(:,chi));
end


% high pass filter lfp_stn and lfp_gp
for chi =1:size(lfp_stn, 2)
    % filter
    lfp_stn(:,chi) = filtfilt(bhp,ahp, lfp_stn(:,chi));
    lfp_gp(:,chi) = filtfilt(bhp,ahp, lfp_gp(:,chi));
end

end



function to_remove_ma = extract_to_remove_ma(idx_noMoveSegs, n)
%% extract the to_remove_ma with MA data

noMovData= zeros(1,n);
for i = 1: size(idx_noMoveSegs, 1)
    noMovData(idx_noMoveSegs(i,1):idx_noMoveSegs(i, 2)) = 1;
end
to_remove_ma = ~noMovData;

end



function power_state_mintime = extract_power_state_mintime(lfp, thr_power, min_samples, fs, freq_band)
%extract_power_state_mintime 
%       get min power state vector by first extracting envelope, setting each
%       power_state to be 1 (lfp_env <= thr_power) or 0 (lfp_env >
%       thr_power1), and get the power_state_mintime based on on
%       power_state and min_samples
%       
%
%       Args:
%           lfp: the lfp vector, 1 * ntemp
%           thr_power: the threshold, a scalar
%           min_samples: an integer scalar
%           fs: sample rate, used for extracting envelope
%           freq_band: band pass filter range, used for extracting envelope
%
%       Return:
%           state_vec: 1 * ntemp vector


% extracting lfp_env through band pass, abs(hilbert) and low pass
lfp_env = lfp;
[blf1,alf1] = butter(2 , 2*freq_band/fs );
lfp_env = filtfilt(blf1,alf1, lfp_env);
lfp_env = abs ( hilbert(lfp_env) ); % hilbert
[blf2,alf2] = butter(2 , 2*0.15/fs );
lfp_env = filtfilt(blf2,alf2, lfp_env); % low pass f


% set power_state =0 for lfp_env > thr_power , 1 for lfp_env <= thr_power
power_state = ones(1,length(lfp_env))*0.5;
power_state(lfp_env> thr_power ) =0;
power_state(lfp_env <= thr_power) =1;

[power_state_mintime] = get_state_mintime( power_state, min_samples );

end




function eye_state_mintime = extract_eye_state_mintime(eye_area_ds, thr_eye_open, time_us, min_samples)
%% extract eye_state_mintime for each time point in based on eye_area_ds and thr_eye_open


% time point for eye and lfp seperately
time_ds = [0:length(eye_area_ds)-1]/15;

% get eye_area with the same time resolution of lfp data using 1-D interpolation
eye_area = interp1(time_ds,eye_area_ds, time_us);


% set eye_state =0 for eye_area < thr_eye_open, 1 for eye_area >= thr_eye_open
eye_state = ones(1,length(eye_area));
eye_state(eye_area < thr_eye_open ) =0;
eye_state(eye_area>= thr_eye_open) =1;

[eye_state_mintime] = get_state_mintime( eye_state, min_samples );

end



function [power_state_mintime] = get_state_mintime( state_in, min_samples )
%GET_STATE_MINTIME 
%       get min state vector based on state_in and min_samples.
%       state_vec is set to be 1/0 if state_in is 1 and the state 1/0 last longer than min_samples, others -1
%                  
%       for example: state_in = [1 1 1 1 0 0 0 1 1], min_samples = 3
%                    =>   state_vec = [1 1 1 1 0 0 0 -1 -1]
%
%       Args:
%           state_in: 1 * ntemp vector
%           
%           min_samples: an integer scalar
%
%       Return:
%           state_vec: 1 * ntemp vector


%%
power_state_mintime = -1*ones(size(state_in));
idx_start = 1; % start index of one state segment
for k = 2:length(state_in)
    
    % state change to a different state (change point) or the last state 
    if ( state_in(k) ~=  state_in(k-1) || k==length(state_in) )
        
        % end index of the last state segment
        idx_end = k-1;
        
        % the last state segement longer than min_elapsed_samples
        if idx_end-idx_start +1 >= min_samples
            
            % previous state is 1, state_vec for this state is assigned 1
            if state_in(idx_start) ==1
                power_state_mintime(idx_start:k-1) = ones(1, k - idx_start )  ;
            end
            
            % previous state is 0, state_vec for this state is assigned 0
            if state_in(idx_start) ==0
                power_state_mintime(idx_start:k-1) = zeros(1, k - idx_start );
            end
        end
        
        idx_start = k;
    end
end

end


function [ vec_segsIndex ] = get_segIndex( states,  min_samples )
% get_segIndex 
%       get segment indices based on state and min_samples. Only get the
%       segment with state == 1 and duration is longer or equall to
%       min_samples
%                  
%       for example: state_in = [1 1 1 1 0 0 0 1 1], min_samples = 2
%                    =>   vec_segsIndex = [1 4; 8 9]
%
%       Args:
%           states: 1 * ntemp vector
%           
%           min_samples: an integer scalar
%
%       Return:
%           vec_segsIndex: nSegs * 2 (idx_start, idx_end)



%%


% find the first index at where state == 1
i = 1;
while(i <=length(states) &&states(i) == 0 )
    i =i+1;
end
idx_state1Start = i;

% extract all the seg start and end indices where state ==1 and duration >= min_samples
vec_segsIndex = [];
for i= idx_state1Start + 1: length(states)
    
    % transition from state 1 to state 0
    if states(i) == 0 && states(i-1) ==1
        idx_state1End = i -1;
        
        % the state 1 segment duration > min_samples
        if idx_state1End - idx_state1Start >= min_samples
            vec_segsIndex = cat(1, vec_segsIndex, [idx_state1Start, idx_state1End]);
        end
    end
    
    % transition from state 0 to state 1
    if states(i) == 1 && states(i-1) == 0
        idx_state1Start = i;
    end
    
    % the last state an its state is 1
    if i == length(states) && states(i) == 1
        idx_state1End = i;
        
        % the state 1 segment duration > min_samples
        if idx_state1End - idx_state1Start >= min_samples
            vec_segsIndex = cat(1, vec_segsIndex, [idx_state1Start, idx_state1End]);
        end
    end
end
end

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

