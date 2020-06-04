clear 
inputfolder = '';

tmin_trial = 1.25;
tmin_reach = 0.5;
tmin_return = 0.5;
t_trial = [0 tmin_trial];
freqs = [6 50];


%% mild
files = dir(fullfile('', '*mild*.mat'));
nfiles = length(files);

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
lfp_trials_mild = lfp_trials;
clear lfp_trials

%% normal extraction
files = dir(fullfile('', '*normal*.mat'));
nfiles = length(files);

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
lfp_trials_normal = lfp_trials;
clear lfp_trials

%% mean lfp_trials: nchns * ntemp * ntrials
lfp_mild = squeeze(mean(mean(lfp_trials_mild, 3),1));
[f_show_mild, power_show_mild ]= fft_show(lfp_mild, fs, freqs);
power_show_mild = norm0to1(power_show_mild);


lfp_normal = squeeze(mean(mean(lfp_trials_normal, 3),1));
[f_show_normal, power_show_normal ]= fft_show(lfp_normal, fs, freqs);
power_show_normal = norm0to1(power_show_normal);


%% plot
linewidth = 3;

color_normal_line = 0.6 * [1 1 1];

color_mild_line = [255,192,203]/255;


figure
hold on
plot(f_show_mild, power_show_mild);
plot(f_show_normal, power_show_normal, 'r');
l_normal = plot(f_show_normal,power_show_normal, 'Color', color_normal_line, 'LineWidth', linewidth);
l_mild = plot(f_show_mild,power_show_mild, 'Color', color_mild_line, 'LineWidth', linewidth);

legend([l_normal, l_mild],'normal', 'mild')
xlabel('frequency')
title('M1 spectrum')

saveas(gcf, 'M1spectrum.png','png')
