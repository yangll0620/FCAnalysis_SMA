function [lfptrials, fs_lfp, T_chnsarea, idxGroups, idxGroupNames] = lfp_goodTrials_align2_J(files, align2, tdur_trial, t_min_reach, varargin)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
%         Args:
%             align2: the event to be aligned (e.g Event.REACHONSET or uint 1-5)
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%             t_minmax_reach, t_minmax_return : min and max reach/return (s) for selecting trials (e.g [0.5 1])
%   
%           Name-Value: 
%               'codesavefolder' - code saved folder
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


coli_align2 = uint32(align2);

coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);



% load fs_lfp and T_chnsarea
file1 = fullfile(files(1).folder, files(1).name);
listOfVariables = who('-file', file1);
if ismember('fs_lfp', listOfVariables)
    load(file1,  'fs_lfp');
else
    load(file1,  'fs');
    fs_lfp = fs;
    clear fs
end
load(file1,  'T_chnsarea');
clear file1 listOfVariables

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    filename = files(filei).name;
    file = fullfile(files(filei).folder, filename);
    
    % load T_idxevent_lfp or T_idxevent
    listOfVariables = who('-file', file);
    if ismember('T_idxevent_lfp', listOfVariables)
        load(file,  'T_idxevent_lfp');
    else
        if ismember('T_idxevent', listOfVariables)
            load(file, 'T_idxevent');
            T_idxevent_lfp = T_idxevent;
            clear T_idxevent
        else
            disp([filename ' not contain T_idxevent_lfp or T_idxevent']);
            continue;
        end
    end
    clear listOfVariables
    
    load(file, 'lfpdata', 'goodTrials');
    
    if filei == 1 % load/extract idxGroups and idxGroupNames once
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
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_min_reach 
            clear t_reach
            continue
        end

               
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
        
        idxdur(1) = idxdur(1) + 1;
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
       
        
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach t_return idxdur lfp_phase_1trial
    end
    
    clear filename file 
    clear('lfpdata', 'T_idxevent_lfp', 'goodTrials');
end
end