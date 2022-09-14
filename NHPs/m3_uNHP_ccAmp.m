function m3_uNHP_ccAmp(animal, varargin)
%
%   Example Usage:
%           m3_uNHP_ccAmp('Kitty', 'newRun', true);
%           m3_uNHP_ccAmp('Jo', 'ei_str', 1, 'ci_str', 1)
%
%   Input:
%       Name-Value: 
%           animal
%           ei_str - event start index
%           ei_end - event start index
%           ci_str - condition start index
%           ci_end - condition end index
%           newRun - true or false(default), for running new or not

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
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
newRun = p.Results.newRun;

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


ccAmpfile_prefix =[animal '-ccAmp'];

%%  input setup

% input folder: extracted raw rest data with grayMatter
if strcmpi(animal, 'Kitty')
    inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');
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


chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
if(strcmpi(animal, 'Jo'))
    cond_cell = [cond_cell {'PD'}];
end

for ei = ei_str: ei_end
    event = EventPhases{ei};
    
    subphsavefolder = fullfile(phsubfolder, event);
    if ~exist(subphsavefolder, 'dir')
        mkdir(subphsavefolder);
    end
    
    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};
        
        disp([codefilename ' ' animal '-' event '-' pdcond])
        
        [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI_extract(event, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
        
        % load(and extract) ccAmpfile
        ccAmpfile = fullfile(savefolder, [ccAmpfile_prefix  '_' pdcond '_' event '.mat']);
        if(~exist(ccAmpfile, 'file') || newRun)
            
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                clear files
                continue;
            end
            
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            if strcmpi(animal, 'Kitty')
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, t_AOI, 't_min_reach', t_minmax_reach(1));
            end
            if strcmpi(animal, 'Jo')
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, t_AOI, 't_min_reach', t_minmax_reach(1)); % lfptrials: nchns * ntemp * ntrials
            end
                        
            % remove unused chns
            ChnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
            lfptrials = lfptrials(ChnsOfI_mask, :, :);
            T_chnsarea = T_chnsarea(ChnsOfI_mask, :);
                        
            %  extract and save ccAmp
            [ccAmp, f_selected] = ccAmp_SKTAllchns_FFT(lfptrials, fs, f_AOI, 'codesavefolder', savecodefolder);
            save(ccAmpfile, 'ccAmp', 'T_chnsarea', 'lfptrials', 'f_selected', 'fs', 'f_AOI');
            
            clear files t_minmax_reach
            clear ChnsOfI_mask
            clear('ccAmp', 'T_chnsarea', 'lfptrials', 'f_selected', 'fs');
        end
    end
end


function [lfptrials, fs_lfp, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, tdur_trial, varargin)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
% Args:
%       files: the files used for extraction       
%       align2: the event to be aligned 
%       tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])     
%
%             
%
%       Name-Value: 
%            't_min_reach' - the minimal reach time used, default 0.5 s
%            't_max_reach' - the max reach time used, default inf
% 
% Return:
%             lfptrials: nchns * ntemp * ntrials
%             chnAreas:
%             fs: sampling rate


% parse params
p = inputParser;
addParameter(p, 't_min_reach', 0.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 't_max_reach', inf, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
t_min_reach = p.Results.t_min_reach;
t_max_reach = p.Results.t_max_reach;


% code start here
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
        if t_reach < t_min_reach || t_reach > t_max_reach
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
            
            if t_reachonset2peakV < 0.3 || t_peakV2reach < 0.3
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


function [lfptrials, fs_lfp, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, tdur_trial, varargin)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately
% [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2PeakV(files, [t_AOI(1) t_AOI(2)], 'codesavefolder', savecodefolder);
%
%   not include trials with t_reach <0.2s
% 
%         Args:
%             align2: the event to be aligned 
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])     
%
%             
%
%       Name-Value: 
%            't_min_reach' - the minimal reach time used, default 0.5 s
%
%   
%           't_max_reach' - the max reach time used, default inf
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


% parse params
p = inputParser;
addParameter(p, 't_min_reach', 0.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 't_max_reach', inf, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
t_min_reach = p.Results.t_min_reach;
t_max_reach = p.Results.t_max_reach;


% code start here
coli_align2 = uint32(align2);
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
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_min_reach  || t_reach > t_max_reach
            clear t_reach
            continue
        end
        
        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
            idx_peakV_ma = idx + idx_reachonset_ma -1;
            t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
            
            if t_reachonset2peakV < 0.3 || t_peakV2reach < 0.3
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
        
        
        % extract phase for 1 trial
        lfp_1trial = lfpdata{tri};
        idxdur = round(tdur_trial * fs_lfp) + idx_time0;
        if idxdur(1) == 0
            idxdur(1) = 1;
        else
            idxdur(1) = idxdur(1) + 1;
        end
        lfp_phase_1trial = lfp_1trial(:,idxdur(1) :idxdur(2));
           
        % cat into lfptrials
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach idxdur lfp_phase_1trial lfp_1trial
    end
end