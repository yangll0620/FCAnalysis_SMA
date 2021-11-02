function [lfptrials, fs_lfp, T_chnsarea, idxGroups, idxGroupNames] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, t_minmax_reach)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately

% 
%         Args:
%             align2: the event to be aligned (e.g Event.REACHONSET or uint 1-5)
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%             t_minmax_reach : min and max reach/return (s) for selecting trials (e.g [0.5 1])
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


coli_align2 = uint32(align2);

coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);
coli_returnonset = uint32(SKTEvent.ReturnOnset);
coli_mouth = uint32(SKTEvent.Mouth);

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'selectedTrials');
    if filei == 1
        load(fullfile(files(filei).folder, filename), 'idxGroups', 'tbl_selectedTrialsMarks');
        idxGroupNames = tbl_selectedTrialsMarks.Properties.VariableNames;
        clear tbl_goodTrialsMarks
    end
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    ntrials = length(lfpdata);
    for tri = 1: ntrials
        
        % ignore trials marked with 0
        if ~selectedTrials(tri)
            continue
        end
        
        % select trials based on reach and return duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_minmax_reach(1) 
            clear t_reach
            continue
        end
        
        
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
        lfp_1trial = lfpdata{tri};
        if idxdur(1) == 0
            idxdur(1) = 1;
        end
        lfp_phase_1trial = lfp_1trial(:, idxdur(1) + 1:idxdur(2));
           
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach idxdur lfp_phase_1trial lfp_1trial
    end
end
end