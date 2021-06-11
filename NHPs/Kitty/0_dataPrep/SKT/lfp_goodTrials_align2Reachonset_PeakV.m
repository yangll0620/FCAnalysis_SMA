function [lfptrials, idxs_peakV, fs_lfp, T_chnsarea, idxGroups, idxGroupNames] = lfp_goodTrials_align2Reachonset_PeakV(files, tdur_trial, t_minmax_reach, t_minmax_return)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
%         Args:
%             align2: the event to be aligned (e.g Event.REACHONSET or uint 1-5)
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%             t_minmax_reach, t_minmax_return : min and max reach/return (s) for selecting trials (e.g [0.5 1])
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
%               
%             idxs_peakV: ntrials * 1, the index point for Peak V
% 
%             chnAreas:
% 
%             fs:
%
%


if nargin < 5
    t_minmax_return = [0 inf];
end
if nargin < 4
    t_minmax_reach = [0 inf];
end


align2 = SKTEvent.ReachOnset;

coli_align2 = uint32(align2);



coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);
coli_returnonset = uint32(SKTEvent.ReturnOnset);
coli_mouth = uint32(SKTEvent.Mouth);

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
idxs_peakV = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'goodTrials', ...
        'smoothWspeed_trial', 'fs_ma', 'T_idxevent_ma');
    if filei == 1
        load(fullfile(files(filei).folder, filename), 'idxGroups', 'tbl_goodTrialsMarks');
        idxGroupNames = tbl_goodTrialsMarks.Properties.VariableNames;
        clear tbl_goodTrialsMarks
    end
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    
    ntrials = size(lfpdata, 3);
    for tri = 1: ntrials
        
        % ignore trials marked with 0
        if ~goodTrials(tri)
            continue
        end       
        
        % select trials based on reach and return duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        t_return = (T_idxevent_lfp{tri, coli_mouth} - T_idxevent_lfp{tri, coli_returnonset}) / fs_lfp;
        if t_reach < t_minmax_reach(1) || t_reach > t_minmax_reach(2)
            clear t_reach
            continue
        end
        if t_return < t_minmax_return(1) || t_reach > t_minmax_return(2)
            clear t_return
            continue
        end
        
      
        
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
           
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        
        % extract the peakV index in the lfptrials
        idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
        idx_reach_ma = T_idxevent_ma{tri, coli_reach};      
        [~, locs] = findpeaks(smoothWspeed_trial(idx_reachonset_ma:idx_reach_ma, tri), 'NPeaks', 1);
        idx_peakV = round((locs / fs_ma - tdur_trial(1) ) * fs_lfp);
        idxs_peakV = [idxs_peakV; idx_peakV];
        clear idx_reachonset_ma idx_reach_ma idx_peakV
        
        
        clear t_reach t_return idxdur lfp_phase_1trial
    end
end
end