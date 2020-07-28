function [lfps_phase_allfiles, lfps_base_allfiles, fs] = lfpm1_allfile(files, eventIdx_phase, tmin_phase, tmax_phase, tphase_bef, tbase_bef)
%% extract the lfp respect to onset of event, [eventonset - tbef eventonset + tmin]
%          and the baseline lfp respect to target onse [targetonset - tbase_bef, targetonse]
%
%
%   return:
%       lfps_phase_allfiles: ntemp * ntrials(ntrials/file * nfiles), average across channels in M1
%       lfps_base_allfiles: ntemp * ntrials(ntrials/file * nfiles), average across channels in M1


depth_M1 = [1.25 1.75] * 8;

lfps_phase_allfiles = [];
lfps_base_allfiles = [];
nfiles = length(files);
for filei = 1 : nfiles
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'fs', 'T_idxevent', 'T_chnsarea');
    
    
    T_idxevent{:,:} = round(T_idxevent{:,:});
    event = T_idxevent{:,eventIdx_phase};
    
    idxmin = tmin_phase * fs;
    idxmax = tmax_phase * fs;
    idxbef = tphase_bef * fs;
    nbase_bef = tbase_bef * fs;
    
    
    
    % extract the date of the exp
    idx = strfind(filename, '_bktdt');
    dateofexp = datetime(filename(idx-6: idx-1), 'InputFormat', 'MMddyy');
    
    % lfp data in m1
    chanUsed_dateofexp = useful_chan_extract('M1', depth_M1, dateofexp);
    chans_m1 = ismember(T_chnsarea.recordingchn, chanUsed_dateofexp);
    lfpm1 = lfpdata(chans_m1, :,:);
    
    lfps_phase_1file = [];
    lfps_base_1file = [];
    for triali = 1 : size(event, 1)
        
        if event(triali, 2) - event(triali, 1) < idxmin || event(triali, 2) - event(triali, 1) > idxmax
            continue;
        end
        
        % phase lfp: nchns * nphasetemp
        lfp_phase = lfpm1(:, event(triali,1)-idxbef : event(triali, 1) + idxmin, triali);
        
        % baseline lfp: nchns * nbasetemp
        lfp_base = lfpm1(:, T_idxevent{triali,1}-nbase_bef : T_idxevent{triali,1}, triali);
        
        
        % lfps_phase_1file: nchns * nphasetemp * ntrials
        lfps_phase_1file = cat(3, lfps_phase_1file, lfp_phase);
        % lfps_base_1file: nchns * nbasetemp * ntrials
        lfps_base_1file = cat(3, lfps_base_1file, lfp_base);
        
        clear lfp_phase lfp_base
    end
    
    
    % if no trial
    if isempty(lfps_phase_1file) && isempty(lfps_base_1file)
        continue;
    end
    
    % lfp_phase_1file (nphasetemp * ntrials/file): average across all the channels 
    lfp_phase_1file = squeeze(mean(lfps_phase_1file, 1));
    lfp_base_1file = squeeze(mean(lfps_base_1file, 1));
    
    % if only one trial
    if(size(lfps_phase_1file, 3) == 1 && size(lfps_base_1file, 3))
        lfp_phase_1file = reshape(lfp_phase_1file, size(lfp_phase_1file, 1) * size(lfp_phase_1file, 2) ,1);
        lfp_base_1file = reshape(lfp_base_1file, size(lfp_base_1file, 1) * size(lfp_base_1file, 2) ,1);
    end
    
    % cat along the dim ntrial
    lfps_phase_allfiles = cat(2, lfps_phase_allfiles, lfp_phase_1file);
    lfps_base_allfiles = cat(2, lfps_base_allfiles, lfp_base_1file);
    
    
    clear lfpdata T_idxevent chans_m1
    clear event idxmin idxmax idxbef lfpm1 
    clear lfp_phase_1file lfp_base_1file
end
