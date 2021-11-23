function [lfptrials, fs, T_chnsarea] = lfp_align2_reachonset(files, tdur_trial, tmin_reach, tmax_reach)

coli_reachonset = 2;
coli_reach = 3;

load(fullfile(files(1).folder, files(1).name),  'fs', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent');
    
    if(height(T_idxevent) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    
    ntrials = size(lfpdata, 3);
    for tri = 1: ntrials
        
        t_reach = (T_idxevent{tri, coli_reach} - T_idxevent{tri, coli_reachonset}) / fs;
        if t_reach < tmin_reach || t_reach > tmax_reach
            continue
        end
            
        idxdur = round(tdur_trial * fs) + T_idxevent{tri, coli_reach};
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
           
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
    end
end