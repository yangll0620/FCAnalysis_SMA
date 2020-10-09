function m2_STKData_narrowfiltered9_11()
%% narrow filtered STK data recorded  in frequency [9 11]Hz
%


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
% band pass frequency
frebp = [9 11];
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');


file_staSTK = fullfile(codecorresParentfolder, 'Jo_SKTTime.mat');
thed_low = 't_10';
thed_high = 't_90';


%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['narrowFiltered' num2str(frebp(1)) '_' num2str(frebp(2)) '-'];

%% starting: narrow filter the lfp data of all the files


%-- stk statistical time --%
load(file_staSTK, 'tbl_staReachT',   'tbl_staReturnT');

files = dir(fullfile(inputfolder, '*.mat'));


% extract brainareas
load(fullfile(files(1).folder, files(1).name), 'T_chnsarea', 'fs');
brainareas = unique(T_chnsarea.brainarea);
nchns_allareas = struct();
for braini = 1: length(brainareas)
    brainarea = brainareas{braini};
    
    
    eval(['nchns_allareas.' brainarea ' = strcmp(T_chnsarea.brainarea, brainarea);'])
end


nfiles = length(files);
f = waitbar(0, ['Narrow Filtering....']);


for filei = 1 : nfiles
    % wait bar
    waitbar(filei/nfiles,f,['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    

    
    % extracted trials in range [-1 tmin_reach+tmin_return + 0.1]s respected to target onset
    tdur_trial = [-1 trange_reach(1) + trange_return(1) + 0.1];
  
    
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
    
        
        
        if idx_touch - idx_reachonset >= trange_reach(1) * fs && idx_touch - idx_reachonset <= trange_reach(2) * fs ...
                && idx_mouth - idx_return >= trange_return(1) * fs && idx_mouth - idx_return <= trange_return(2) * fs ...
                && idx_mouth - idx_target <= tmax_trial * fs &&  idx_mouth - idx_target >= tdur_trial(2) * fs % satisfied trials: in trange_reach and trange_return and less than tdur_trial
            
            
        end
        
        
        idxs_trial = round(tdur_trial * fs) + idx_target;
        if idxs_trial(1) == 0
            idxs_trial(1) = 1;
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
        
        clear ntrial
        
        clear idx_target idx_touch idx_mouth
        
    end
    
    
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
    
    
    %%% band pass filter
    [nchns, ntemps, ntrials] = size(lfpdata);
    filterdlfp = zeros(nchns, ntemps, ntrials);
    for chni = 1 : nchns
        for triali = 1: ntrials
            filterdlfp(chni,:, triali) = filter_bpbutter(lfpdata(chni,:,triali),frebp,fs);
        end
    end
    
    
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
    [nchns, ntemp, ntrials] = size(filterdlfp);
    lfpdata = zeros(ntemp, nchns, ntrials);
    for chi = 1 : nchns
        lfpdata(:, chi, :) = filterdlfp(chi, :, :);
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
close(f);
disp(['narrow filtered lfpdata  are saved to ' savefolder])

