function [lfptrials, fs_lfp, T_chnsarea, idxGroups, idxGroupNames] = lfpseg_selectedTrials_align2PeakV(files, tdur_trial, varargin)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately
% [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2PeakV(files, [t_AOI(1) t_AOI(2)], 'codesavefolder', savecodefolder);
% 
%         Args:
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);


load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'selectedTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
    
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
        
        % find peakV and its timepoint
        idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
        idx_reach_ma = T_idxevent_ma{tri, coli_reach};
        [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
        idx_peakV_ma = idx + idx_reachonset_ma -1;
        t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
        t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;

        if t_reachonset2peakV < abs(tdur_trial(1)) || t_peakV2reach < tdur_trial(2)
            clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma 
            clear t_reachonset2peakV  t_peakV2reach
            continue;
        end
        
        % extract trial with t_dur
        idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
        idxdur = round(tdur_trial * fs_lfp) + idx_peakV_lfp;
        lfp_1trial = lfpdata{tri};
        if idxdur(1) == 0
            idxdur(1) = 1;
        else
            idxdur(1) = idxdur(1) + 1;
        end
        lfp_phase_1trial = lfp_1trial(:,idxdur(1) :idxdur(2));
           
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach idxdur lfp_phase_1trial lfp_1trial
    end
end
end