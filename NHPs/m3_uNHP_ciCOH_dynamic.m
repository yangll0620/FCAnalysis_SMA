function m3_uNHP_ciCOH_dynamic(animal, varargin)
%
%   Example Usage:
%           m3_uNHP_ciCOH_dynamic('Kitty', 'newRun', false, 'shuffleN_psedoTest', 200);
%           m3_uNHP_ciCOH_dynamic('Jo', 'newRun', false, 'shuffleN_psedoTest', 200);
%
%   Input:
%       Name-Value: 
%           animal
%           ei_str - event start index
%           ei_end - event start index
%           ci_str - condition start index
%           ci_end - condition end index
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500

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

cond_cell = {'normal', 'PD'};

% parse params
p = inputParser;
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', length(cond_cell), @isscalar);
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));

parse(p,varargin{:});
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
newRun = p.Results.newRun;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;

% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);

if strcmpi(animal, 'Kitty')
    NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'SKT','fs500Hz', 'longerTrials', codefilename);
end

if strcmpi(animal, 'Jo')
    NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'SKT','fs500Hz', codefilename);
end
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);


%% save setup
savefolder = codecorresfolder;
codesavefolder = fullfile(savefolder, 'code');
if exist(codesavefolder, 'dir')
    rmdir(codesavefolder,'s');
end
copyfile2folder(codefilepath, codesavefolder);


ciCohfile_prefix =[animal '-ciCoh_Dynamic_'];

%%  input setup

% input folder: extracted raw rest data with grayMatter
if strcmpi(animal, 'Kitty')
    inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');
end
if strcmpi(animal, 'Jo')
    inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');
end


f_AOI = [8 40];

align2 = SKTEvent.ReachOnset;
if strcmpi(animal, 'Kitty')
    t_trialdur = [-1.5 1.5];
end
if strcmpi(animal, 'Jo')
    t_trialdur = [-1.2 1.5];
end
t_min_reach = 0.5;

t_win = 0.2;
t_overlap = 0.15;
t_shift = t_win - t_overlap;

%% Code start here

chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', codesavefolder);


for ci = ci_str : ci_end
    pdcond = cond_cell{ci};
    
    disp([codefilename ' ' animal  '-' pdcond])
    
    % load(and extract) ciCohfile
    ciCohfile = fullfile(savefolder, [ciCohfile_prefix  '_' pdcond '.mat']);
    
    if(~exist(ciCohfile, 'file') || newRun)
        
        if strcmp(pdcond, 'PD')
            files1 = dir(fullfile(inputfolder, ['*_mild_*.mat']));
            files2 = dir(fullfile(inputfolder, ['*_moderate_*.mat']));
            files = [files1; files2];
            clear files1 files2
        else
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
        end
        
        
        if isempty(files)
            clear files
            continue;
        end
        
        if strcmpi(animal, 'Kitty')
            [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, t_trialdur, 't_min_reach', t_min_reach);
        end
        if strcmpi(animal, 'Jo')
            [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, t_trialdur, 't_min_reach', t_min_reach); % lfptrials: nchns * ntemp * ntrials
        end
        
        % remove unused chns
        ChnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
        lfptrials = lfptrials(ChnsOfI_mask, :, :);
        T_chnsarea = T_chnsarea(ChnsOfI_mask, :);
        
        %  extract and save ccAmp along time
        ciCoh = [];
        t_selected = [];
        lfptrials_win = []; % lfptrials_win: nchns * ntemp * ntrials * nts
        n_shift = round(t_shift * fs);
        n_win = round(t_win * fs);
        len = size(lfptrials, 2);
        for idx_str =  1 : n_shift : len - n_win + 1
            idx_end = idx_str + n_win - 1;
            
            t = ((idx_str -1 + idx_end)/2)/fs + t_trialdur(1);
            lfptrials_win = cat(4, lfptrials_win, lfptrials(:, idx_str:idx_end, :));

            t_selected = [t_selected; t];
            clear idx_end
        end
        clear idx_str
        
        for ti = 1 : size(lfptrials_win, 4)
            lfp = squeeze(lfptrials_win(:, :, :, ti));
            [ciCoh_phase, f_selected] = ciCohSKTAllchns_FFT_NoAmp(lfp, fs, f_AOI, 'codesavefolder', codesavefolder);
            ciCoh = cat(4, ciCoh, ciCoh_phase);
            clear lfp ciCoh_phase
        end
        
        save(ciCohfile, 'ciCoh', 'T_chnsarea', 'lfptrials', 'f_selected', 't_selected', 'fs', 'f_AOI', 't_trialdur', 'lfptrials_win');
        
        clear files ChnsOfI_mask
        clear('ciCoh', 'T_chnsarea', 'lfptrials', 'f_selected', 'fs', 't_selected', 'lfptrials_win');
    end
    
    %%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
    load(ciCohfile, 'lfptrials_win', 'fs', 'f_AOI');
    for ti = 1 : size(lfptrials_win, 4)
        lfp = squeeze(lfptrials_win(:, :, :, ti));
        psedociCoh_extract_save(shuffleN_psedoTest, lfp, fs, f_AOI, ciCohfile, codesavefolder, ti);
        clear lfp
    end
    clear lfptrials_win fs 

    
    clear pdcond ciCohfile
end



function psedociCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohfile, codesavefolder, ti)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohPhasefile:
%       
%       Name-Value: 
%           'codesavefolder' - code saved folder         
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile


psedoVar = whos('-file',ciCohfile, 'psedociCohs');
if(isempty(psedoVar))
    psedociCohs = struct();
else
    load(ciCohfile, 'psedociCohs');
end
clear psedoVar


fiename = ['t' num2str(ti)];
if ~isfield(psedociCohs, fiename)
    psedociCohs.(fiename) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs.(fiename), 4) + 1;
end

if shuffi_str > suffi_end
    return;
end

disp(['ti = ' num2str(ti) ' psedo test start at ' num2str(shuffi_str) ' times'])
nchns = size(lfptrials, 1);
for si = shuffi_str : suffi_end
    for chni = 1 : nchns -1
        lfp1 = squeeze(lfptrials(chni, :, :));
        for chnj = chni + 1 : nchns
            lfp2 = squeeze(lfptrials(chnj, :, :));
            [psedociCoh, ~] = psedo_ciCohSKT_FFT_NoAmp(lfp1, lfp2, fs, f_AOI, 'codesavefolder', codesavefolder);
            
            if chni == 1 && chnj == chni + 1
                nf = size(psedociCoh, 1);
                psedociCohs1Time = zeros(nchns, nchns, nf, 1);
                clear nf
            end
            
            psedociCohs1Time(chni, chnj, :, :) = psedociCoh;
            clear lfp2 psedociCoh
        end
        clear lfp1
    end
    psedociCohs.(fiename) = cat(4, psedociCohs.(fiename), psedociCohs1Time);
    clear psedociCohs1Time
    if(mod(si,10) == 0)
        disp(['ti = ' num2str(ti) ' psedo test finished ' num2str(si) ' times'])
        save(ciCohfile, 'psedociCohs', '-append');
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
        if(idxdur(1) <= 0)
            continue;
        end
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
addParameter(p, 't_min_reach', 0, @(x) assert(isnumeric(x) && isscalar(x)));
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