function m2_STKData_seg()
%% seg each trial to equal length
%
% different for normal, mild and moderate
%  
% tdur_trial = [-1, tmin_reach+t_min_return+0.1] respect to targetonset
%
%
%   a. average across area if not DBS
%
% Output:
%       ntemp * nchns * ntrials (one/area + bipolar DBS)

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
tmp = char(regexp(codefilepath, 'NHPs/[A-Za-z]*', 'match'));
animal = tmp(length('NHPs/')+1:end);


%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');

% file for statistical STK performance
file_staSTK = fullfile(codecorresParentfolder, 'm1_SKTData_timeStatiscal', [animal '_SKTTime.mat']);



%% save setup
savefolder = codecorresfolder;


%% code start here

pdcond = 'normal';
files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
m2_STKData_seg_perCond(pdcond, file_staSTK, files, savefolder, animal);


pdcond = 'mild';
files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
m2_STKData_seg_perCond(pdcond, file_staSTK, files, savefolder, animal);


function m2_STKData_seg_perCond(pdcond, file_staSTK, files, savefolder, animal)


savefilename_addstr = 'segsamelength';
load(file_staSTK, 'tbl_staReachT',   'tbl_staReturnT');
thed_low = 't_10';
thed_high = 't_90';

tmax_trial = 3;


% extract trange_reach and trange_return
if strcmp(pdcond, 'normal')
    trange_reach = [tbl_staReachT{'normal', thed_low} tbl_staReachT{'normal', thed_high}];
    trange_return = [tbl_staReturnT{'normal', thed_low} tbl_staReturnT{'normal', thed_high}];
else
    if strcmp(pdcond, 'mild')
        trange_reach = [tbl_staReachT{'mild', thed_low} tbl_staReachT{'mild', thed_high}];
        trange_return = [tbl_staReturnT{'mild', thed_low} tbl_staReturnT{'mild', thed_high}];
    else
        if strcmp(pdcond, 'moderate')
            trange_reach = [tbl_staReachT{'moderate', thed_low} tbl_staReachT{'moderate', thed_high}];
            trange_return = [tbl_staReturnT{'moderate', thed_low} tbl_staReturnT{'moderate', thed_high}];
        end
    end
end

% extracted trials in range [-1 tmin_reach+tmin_return + 0.1]s respected to target onset
tdur_trial = [-1 trange_reach(1) + trange_return(1) + 0.1];


% extract brainareas
load(fullfile(files(1).folder, files(1).name), 'T_chnsarea');
brainareas = unique(T_chnsarea.brainarea);
nchns_allareas = struct();
for braini = 1: length(brainareas)
    brainarea = brainareas{braini};
    
    
    eval(['nchns_allareas.' brainarea ' = strcmp(T_chnsarea.brainarea, brainarea);'])
end


%% inside for
nfiles = length(files);
for filei = 1: nfiles

    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    [m, n, l] = size(lfpdata);
    if m==1 || n==1 || l==1
        continue;
    end
    
    disp([num2str(filei) '/' num2str(nfiles) ': ' filename])
    
    
    %%% extract lfp_trials satifiy the critirion
    lfp_trials = struct();
    for braini = 1: length(brainareas)
        brainarea = brainareas{braini};
        
        eval(['lfp_trials.' brainarea ' = [];']);
    end
    for triali = 1: height(T_idxevent)
        
        idx_target = T_idxevent{triali,'TargetTime'};
        idx_reachonset = T_idxevent{triali, 'ReachTimeix'};
        idx_touch = T_idxevent{triali, 'TouchTimeix'};
        idx_return = T_idxevent{triali, 'ReturnTimeix'};
        idx_mouth = T_idxevent{triali, 'MouthTimeix'};
        
        
        if idx_touch - idx_reachonset < trange_reach(1) * fs || idx_touch - idx_reachonset > trange_reach(2) * fs ...
                && idx_mouth - idx_return < trange_return(1) * fs && idx_mouth - idx_return > trange_return(2) * fs ...
                && idx_mouth - idx_target > tmax_trial * fs &&  idx_mouth - idx_target < tdur_trial(2) * fs
            
            continue;
            % not satisfied trials: in trange_reach and trange_return and less than tdur_trial
        end
        
        % indices for extracted trial
        idxs_trial = round(tdur_trial * fs) + idx_target;
        if idxs_trial(1) == 0
            idxs_trial(1) = 1;
        else
            continue;
        end
        
        
        % combine to each brain area field
        for braini = 1: length(brainareas)
            brainarea = brainareas{braini};
            
            
            % extract the chns for brainarea
            eval(['chns = nchns_allareas.' brainarea ';'])
            % extract lfp1trial for brainarea
            lfp1trial = lfpdata(chns, idxs_trial(1): idxs_trial(2), triali);
            
            
            % combined to particular brainarea field of lfp_trails
            eval(['lfp_trials.' brainarea '  = cat(3, lfp_trials.' brainarea ', lfp1trial);'])
            
            clear brainarea lfp1trial chns
        end
        
        
        clear idx_target idx_reachonset idx_touch idx_return idx_mouth
        
    end % end for triali
    
    
    
    %%% average across area if not DBS, extract lfpdata and T_chnsarea_new
    lfpdata = [];
    T_chnsarea_new = T_chnsarea([], :);
    for braini = 1: length(brainareas)
        brainarea = brainareas{braini};
        
        eval(['lfp = lfp_trials.' brainarea ';']);
        
        if ~any(strcmp(brainarea, {'STN', 'GP'}))
            lfp = mean(lfp, 1); 
        end
        
        % concatenate to lfpdata
        lfpdata = cat(1, lfpdata, lfp);
        
        
        %%% new_T_chnsarea
        
        % extract the chns for brainarea
        eval(['chns = nchns_allareas.' brainarea ';'])
        
        T_area = T_chnsarea(chns, :);
        if ~any(strcmp(brainarea, {'STN', 'GP'})) % only has one row if not DBS
            T_area = T_area(1,:);
            T_area.recordingchn = NaN;
        end
        T_chnsarea_new = [T_chnsarea_new; T_area];
        clear T_area
        
        clear lfp brainarea
    end
    T_chnsarea = T_chnsarea_new;
    clear T_chnsarea_new
    
    
    
    % extract chnAreas cell for used in python
    chnAreas =T_chnsarea.brainarea;
    idx_STN = find(strcmp(T_chnsarea.brainarea, 'STN'));
    for i = 1: length(idx_STN)
        chnAreas{idx_STN(i)} = ['stn' num2str(i-1) '-' num2str(i)];
    end
    idx_GP = find(strcmp(T_chnsarea.brainarea, 'GP'));
    for i = 1: length(idx_GP)
        chnAreas{idx_GP(i)} = ['gp' num2str(i-1) '-' num2str(i)];
    end
    
    
    % save
    
    % change to ntemp * nchns * ntrials
    lfpdata_tmp = lfpdata;
    clear lfpdata
    [nchns, ntemp, ntrials] = size(lfpdata_tmp);
    lfpdata = zeros(ntemp, nchns, ntrials);
    for chi = 1 : nchns
        lfpdata(:, chi, :) = lfpdata_tmp(chi, :, :);
    end
    
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ...
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'chnAreas');
    
    
    clear lfpdata fs T_chnsarea T_idxevent
    clear nmin_trial  nmin_reach nmin_return
    clear lfp_trials filterdlfp
    clear chnAreas idx_GP idx_STN
end