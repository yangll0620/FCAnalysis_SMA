clear 
inputfolder = '';
files = dir(fullfile('', '*mild*.mat'));
nfiles = length(files);


tmin_trial = 1.25;
tmin_reach = 0.5;
tmin_return = 0.5;
t_trial = [0 tmin_trial];

% area chns
load('chans_m1.mat', 'chans_m1');
chns = chans_m1;

lfp_trials = [];
for filei = 1 : nfiles
    
    % load data
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    nmin_trial = round(tmin_trial * fs);
    nmin_reach = round(tmin_reach * fs);
    nmin_return = round(tmin_return * fs);
    
    for triali = 1: height(T_idxevent)
        idx_target = T_idxevent{triali,'TargetTime'};
        idx_touch = T_idxevent{triali, 'TouchTimeix'};
        idx_mouth = T_idxevent{triali, 'MouthTimeix'};
        
        
        if idx_mouth - idx_target > nmin_trial && idx_touch - idx_target > nmin_reach && idx_mouth - idx_touch > nmin_return
            n_trial = round(t_trial * fs) + idx_target;
            
            lfp1trial = lfpdata(chns, n_trial(1): n_trial(2), triali);
            lfp_trials  = cat(3, lfp_trials, lfp1trial);
            
            clear ntrial lfp1trial
        end
        
        clear idx_target idx_touch idx_mouth
        
    end
    clear filename lfpdata T_chnsarea T_idxevent
    clear nmin_trial nmin_reach nmin_return triali
end


%% mean lfp_trials: nchns * ntemp * ntrials
lfp = squeeze(mean(mean(lfp_trials, 3),1));
lfp_meanarea = squeeze(mean(lfp_trials, 1));

freqs = [6 50];

%% fft
for triali = 1:size(lfp_meanarea, 2)
    fft_show(lfp_meanarea(:,triali), fs);
end

