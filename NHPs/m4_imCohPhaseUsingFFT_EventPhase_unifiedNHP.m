function m4_imCohPhaseUsingFFT_EventPhase_unifiedNHP(animal, varargin)
%
%   Example Usage:
%           m4_imCohPhaseUsingFFT_EventPhase_unifiedNHP('Jo', 'ei_str', 1, 'ci_str', 1)
%
%   Input:
%       Name-Value: 
%           animal
%           ei_str - event start index
%           ei_end - event start index
%           ci_str - condition start index
%           ci_end - condition end index

codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));

cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract(animal);

% parse params
p = inputParser;
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', length(EventPhases), @isscalar);
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', length(cond_cell), @isscalar);
parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;

% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);

NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'SKT','fs500Hz', codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, 'dir')
    rmdir(savecodefolder,'s');
end
copyfile2folder(codefilepath, savecodefolder);

phsubfolder = fullfile(savefolder, 'phases');
if ~exist(phsubfolder, 'dir')
    mkdir(phsubfolder);
end


ciCohPhasefile_prefix =[animal ' ciCohPhasefile'];

%%  input setup

% input folder: extracted raw rest data with grayMatter
if strcmpi(animal, 'Kitty')
    inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_goodReach');
end
if strcmpi(animal, 'Jo')
    inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');
end


f_AOI = [8 40];


shuffleN_psedoTest = 500;


%% Code start here
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal, 'codesavefolder', savecodefolder);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', savecodefolder);

unwanted_DBS = unwanted_DBS_extract(animal,'codesavefolder', savecodefolder);
noisy_chns = noisy_chns_extract(animal,'codesavefolder', savecodefolder);
notAOI_chns = notInterested_chns_extract(animal, 'codesavefolder', savecodefolder);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns

for ei = ei_str: ei_end
    event = EventPhases{ei};
    
    subphsavefolder = fullfile(phsubfolder, event);
    if ~exist(subphsavefolder, 'dir')
        mkdir(subphsavefolder);
    end
    
    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};
        
        disp([codefilename ' ' animal '-' event '-' pdcond])
        
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, 'codesavefolder', savecodefolder);
        
        % load(and extract) ciCohPhasefile
        ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond '_' event '_align2' align2name '.mat']);
        if(~exist(ciCohPhasefile, 'file'))
            
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end
            
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            if strcmpi(animal, 'Kitty')
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder',savecodefolder);
            end
            if strcmpi(animal, 'Jo')
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder); % lfptrials: nchns * ntemp * ntrials
            end
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            
            %  extract and save deltaphis_allChnsTrials and cicoh
            [deltaphis_allChnsTrials, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfptrials, fs, f_AOI, 'codesavefolder', savecodefolder);
            

            ntrials = size(lfptrials, 3);
            save(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
                      
            %  extract and save psedociCohs
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
            
            
            clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
            clear t_minmax_reach files lfptrials
        end
        load(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
        if ~exist('psedociCohs','var')
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if strcmpi(animal, 'Kitty')
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder);
            end
            if strcmpi(animal, 'Jo')
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder); % lfptrials: nchns * ntemp * ntrials
            end
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
            load(ciCohPhasefile, 'psedociCohs');
            clear t_minmax_reach lfptrials
        end
        
        nshuffle = size(psedociCohs, 4);
        if(nshuffle < shuffleN_psedoTest)
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder);
            else
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder); % lfptrials: nchns * ntemp * ntrials
            end
            
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
            load(ciCohPhasefile, 'psedociCohs');
            clear t_minmax_reach lfptrials
        end
         
        clear pdcond subpdsavefolder align2 t_AOI align2name ciCohPhasefile
        clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
    end
end


function [lfptrials, fs_lfp, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach, varargin)
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

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    filename = files(filei).name;
    file = fullfile(files(filei).folder, filename);
    
    load(file, 'lfpdata', 'T_idxevent_lfp', 'goodTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');

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
        if t_reach < t_minmax_reach(1) 
            clear t_reach
            continue
        end


        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial(idx_reachonset_ma: idx_reach_ma, tri));
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
            idx_time0 = idx_peakV_lfp;
            
            clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
            clear t_reachonset2peakV  t_peakV2reach
            clear idx_peakV_lfp
        else
            idx_time0 = T_idxevent_lfp{tri, coli_align2}; 
        end

               
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs_lfp) + idx_time0;
        idxdur(1) = idxdur(1) + 1;
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
       
        
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach t_return idxdur lfp_phase_1trial idx_time0
    end
    
    clear filename file ntrials
    clear('lfpdata', 'T_idxevent_lfp', 'goodTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
end

