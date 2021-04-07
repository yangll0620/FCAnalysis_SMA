function m1_restData_cleaned_extract()
%% extract cleaned restData using files Ying used
%
%   1. remove data of channels in m1array_to_remove
%   
%   2. remove data with eye-closed
%
%   3. remove data marked with movement 
%   
%   Inputs:
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
[datafolder, ~, ~, ~] = exp_subfolders();

%% animal name
animal = 'Bug';

%% input setup

% input folder: root2 in server
folder_processed_root2 = fullfile('/home/lingling/root2/Animals2', animal, 'Recording/Processed/DataDatabase');


% file_sycChecked by Ying
file_sycChecked = fullfile(datafolder, 'root2', 'Ying Yu', 'Bug data','SycAllFileChecked.mat');

% threshold used by Ying to extract cleaned resting data
configFile = fullfile(datafolder, 'config_m1lf_fromYing.mat');

% master file
file_masterxls = fullfile(datafolder, animal, 'BugMasterDatabase.xlsx'); 


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

% load configFile
load(configFile, 'config_m1lf');


% load file SycAllFileChecked.mat
load(file_sycChecked, 'SycAllFileChecked');
% only use dates marked with SycAllFileChecked.checkstatus == [1, 1]
masks = (SycAllFileChecked.checkstatus(:, 1) == 1 & SycAllFileChecked.checkstatus(:, 2) == 1);
dates_used = unique(SycAllFileChecked.date(masks,:));


% load master file as a table
t_master = readtable(file_masterxls);
% string column OutputFolderName and BriefDescription
t_master.OutputFolderName = string(t_master.OutputFolderName);
t_master.BriefDescription = string(t_master.BriefDescription);


nfiles = length(dates_used);
for i = 1: nfiles
    
    idx = find(t_master.OutputFolderName == dates_used{i} & t_master.BriefDescription == 'Resting');
    
    if (length(idx) ~= 1) % idx = [] or idx has more than 1 value
        
        disp(['Date: '   dates_used{i} ' has ' num2str(length(idx)) ' file, skipped.'])
        continue;
    end
    
    
    % tdt block for resting on date_used
    tdtblock = t_master.TDTBlock(idx);
    
    % sync file pattern
    folder_allsync = fullfile(folder_processed_root2, dates_used{i}, ['Block-' num2str(tdtblock)]);
    filepattern_allsync = fullfile(folder_allsync, [animal '*_all_sync.mat']);
    files_allsync = dir(filepattern_allsync);
    
    % all sync file not exist or there are more than one sync file, skip
    if(length(files_allsync) ~= 1)
        disp([filepattern_allsync ' has ' num2str(length(files_allsync)) ' files, not 1 file, skipped.'])
        continue;
    end
    
    % the all sync file 
    file_allsync = fullfile(files_allsync.folder, files_allsync.name);
    
    
    %%%  extract the segments %%%
    % get the thr_power
    date_num = datenum(strrep(dates_used{i}, [animal '_'], ''), 'mmddyy');
    condition = parsePDCondition(date_num, animal);

    
    if strcmp(condition, 'mild') || strcmp(condition, 'moderate')
        eval(['thr_power = config_m1lf.max_' lower(animal) '_pd;'])
    else
        if strcmp(condition, 'normal')
            eval(['thr_power = config_m1lf.max_' lower(animal) '_normal;'])
        else
            disp([dates_used{i} ', Ignore condition = ' condition])
            
            continue;
        end
    end
    
    disp([num2str(i) '/' num2str(nfiles) ': extracting segments from ' dates_used{i} '-Block' num2str(tdtblock)])
    
    %
    [data_segments, segsIndex, therapy, fs] = segments_extract_fromAllSyncMat(file_allsync, arraychns_to_remove, thr_power, thr_eye_open, t0, min_comp_time, config_m1lf.freq_band);
    
    
    if(isempty(data_segments))
        continue;
    end
    
    str_therapy = '';
    if (strcmp(condition, 'mild') || strcmp(condition, 'moderate')) && ~strcmp(therapy, 'off')
            str_therapy = therapy;
    end
    
    savefile = fullfile(savefolder, [savefilename_prefix  condition str_therapy '_'  datestr(date_num,'yyyymmdd') '_tdt'  num2str(tdtblock) '.mat']);
    save(savefile, 'data_segments', 'segsIndex', 'fs')
    
    clear allsyncfile dateblockstr thr_power1
    clear tmp condition 
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
lfp_array = data.lfp_array;
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

