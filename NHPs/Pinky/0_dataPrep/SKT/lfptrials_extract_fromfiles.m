function lfp_trials = lfptrials_extract_fromfiles(files, t_trial, tmin_trial, tmin_reach, tmin_return)
%% extract the combined lfp trials from all the files 
%
%       input:
%           files = dir(fullfile(inputfolder, ['*_' pdcond '*.mat']));
%
%           t_trial: a vector with two time values respect to target on set e.g[-1 1]
%
%           tmin_trial, tmin_reach, tmin_return: min t_trial, t_reach t_return values; 
%                                                ignore the trial less than any of the condition
%
%       return:
%           lfp_trials: extracted all trials struct (e.g lfp_trials.M1: nchns * ntemp * ntrials)


% extract brain areas from files(1)
load(fullfile(files(1).folder, files(1).name), 'T_chnsarea', 'fs');
brainareas = unique(T_chnsarea.brainarea);


% init lfp_trails struct
lfp_trials = struct();
lfp_trials.fs = fs;
nchns_allareas = struct();
for braini = 1: length(brainareas)
    brainarea = brainareas{braini};
    
    eval(['lfp_trials.' brainarea ' = [];']);
    
    eval(['nchns_allareas.' brainarea ' = strcmp(T_chnsarea.brainarea, brainarea);'])
end


nfiles = length(files);
for filei = 1 : nfiles
    
    % load data
    load(fullfile(files(filei).folder, files(filei).name), 'lfpdata', 'fs', 'T_idxevent');
    
    
    nmin_trial = round(tmin_trial * fs);
    nmin_reach = round(tmin_reach * fs);
    nmin_return = round(tmin_return * fs);
    
    
    for triali = 1: height(T_idxevent)
        idx_target = T_idxevent{triali,'TargetTime'};
        idx_touch = T_idxevent{triali, 'TouchTimeix'};
        idx_mouth = T_idxevent{triali, 'MouthTimeix'};
        
        
        if idx_mouth - idx_target > nmin_trial && idx_touch - idx_target > nmin_reach && idx_mouth - idx_touch > nmin_return
            
            n_trial = round(t_trial * fs) + idx_target;
            if n_trial(1) < 1
                n_trial(1) = 1;
            end
            
            % combine to each brain area field
            for braini = 1: length(brainareas)
                brainarea = brainareas{braini};
                
                
                % extract the chns for brainarea
                eval(['chns = nchns_allareas.' brainarea ';'])
                % extract lfp1trial for brainarea
                lfp1trial = lfpdata(chns, n_trial(1): n_trial(2), triali);
                
                
                % combined to particular brainarea field of lfp_trails 
                eval(['lfp_trials.' brainarea '  = cat(3, lfp_trials.' brainarea ', lfp1trial);'])
                
                clear brainarea lfp1trial chns
            end
            
            clear ntrial 
        end
        
        clear idx_target idx_touch idx_mouth
        
    end
    
    clear lfpdata fs T_idxevent
    clear nmin_trial nmin_reach nmin_return
    clear triali
end