function m4_fs500Hz_FreezeSegs_extract()
% extract freeze segments
%
%


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'toolbox')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);


%%  input setup
inputfolder_FreezeEpisode = fullfile(codecorresParentfolder, 'm3_fs500Hz_freezeSKTData_EpisodeExtract');

seg_tseg = 0.2;
pdcond = 'moderate';


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

savefile_prefix = [animal '-FreezeSegs'];


%% Code start here
files_freezeEpisode = dir(fullfile(inputfolder_FreezeEpisode, '*moderate*.mat'));



% lfpsegments of freeze phases using only reach freeze trials
[lfpsegs_Freeze, fs, T_chnsarea]= seg2Short_FreezeTrials_Segments(files_freezeEpisode, 'tseg', seg_tseg);


% lfpsegments of reach trials
ePhases = {'preMove'; 'earlyReach';  'PeakV'; 'lateReach'};
inputfolder_ReachTrials = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');
files_ReachTrials = dir(fullfile(inputfolder_ReachTrials, ['*' pdcond '_*.mat']));
for ei = 1 : length(ePhases)
    ePhase = ePhases{ei};
    [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI_extract(ePhase, animal, pdcond, 'codesavefolder', savecodefolder);
    [lfptrials, ~, ~] = lfptrials_animalK_align2(files_ReachTrials, align2, t_AOI, ...
            't_min_reach', 0.5, 't_max_reach', inf);

    lfptrials_Reach.(ePhase) = lfptrials;

    clear ePhase align2 t_AOI lfptrials
end


savefilename = [savefile_prefix '-ReachTrials-FreezeSegs.mat'];
savefile = fullfile(savefolder, savefilename);
save(savefile, 'lfpsegs_Freeze', 'lfptrials_Reach', 'fs', 'T_chnsarea', 'seg_tseg');
clear lfpsegs_freeze fs T_chnsarea savefilename savefile


end


function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_animalK_align2(files, align2, tdur_trial, varargin)
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
end


function [lfpsegs_freeze, fs, T_chnsarea]= seg2Short_FreezeTrials_Segments(files, varargin)
%
%  using Freeze trials to extract lfp segments
%
%   Input:
%       files: 
%
%       freezeType: 'InitFreeze', 'ReachFreeze', or 'ManiFreeze'
%    
%       Name-Value: 
%           
%           tseg: the seg duration in second, default = 0.2
%
%   Return:
%         lfpsegs_freeze = 
%         
%           struct with fields:
%         
%               earlyFreeze: [3×100×380 double]
%              middleFreeze: [3×100×4326 double]
%                lateFreeze: [3×100×380 double]
%             afterFreeze200ms: [3×100×380 double]
%            
%       fs, T_chnsarea


if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpsegs_freeze = [];
    fs = [];
    T_chnsarea = [];

    return;
end

% parse params
p = inputParser;
addParameter(p, 'tseg', 0.2, @isscalar);

parse(p,varargin{:});
tseg = p.Results.tseg;


lfpsegs_freeze = struct();
lfpsegs_freeze.InitFreeze.earlyFreeze = [];
lfpsegs_freeze.InitFreeze.middleFreeze = [];
lfpsegs_freeze.InitFreeze.lateFreeze = [];
lfpsegs_freeze.InitFreeze.afterFreeze200ms = [];
lfpsegs_freeze.ReachFreeze.earlyFreeze = [];
lfpsegs_freeze.ReachFreeze.middleFreeze = [];
lfpsegs_freeze.ReachFreeze.lateFreeze = [];
lfpsegs_freeze.ReachFreeze.afterFreeze200ms = [];
lfpsegs_freeze.ManiFreeze.earlyFreeze = [];
lfpsegs_freeze.ManiFreeze.middleFreeze = [];
lfpsegs_freeze.ManiFreeze.lateFreeze = [];
lfpsegs_freeze.ManiFreeze.afterFreeze200ms = [];



for fi = 1: length(files)
    loadfilename = files(fi).name;
    load(fullfile(files(fi).folder, loadfilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'freezStruct', 'selectedTrials');
    
    % check fs and T_chnsarea same across files
    if ~exist('fs_unit', 'var')
        fs_unit = fs_lfp;
    else
        if(fs_unit ~= fs_lfp)
            disp(['fs ~= fs_lfp for file ' loadfilename])
            continue;
        end
    end
    
    if ~exist('T_chnsarea_unit', 'var')
        T_chnsarea_unit = T_chnsarea;
    else
        if(~isequal(T_chnsarea_unit, T_chnsarea))
            disp(['T_chnsarea_unit ~= T_chnsarea for file ' loadfilename])
        end
    end
    
    nseg = round(tseg * fs_unit);
    
    freezEpisodes = freezStruct.freezEpisodes;
    for frzi = 1 : length(freezEpisodes) 
        
        %%%  extract freeze lfp data: seglfp
        tri = freezEpisodes{frzi}.triali;
        if ~selectedTrials(tri)
            clear tri
            continue;
        end

        freezeUsedType = freezEpisodes{frzi}.freezeType;
        if contains(freezeUsedType, 'init Move') 
            freezeType = 'InitFreeze';
        elseif contains(freezeUsedType, 'Reach')
            freezeType = 'ReachFreeze';
        elseif contains(freezeUsedType, 'Manipulation')
            freezeType = 'ManiFreeze';
        end
        clear freezeUsedType
        
        
        % Early Freeze
        t_start_Early = freezEpisodes{frzi}.freezeTPhaseS(1);
        if strcmpi(freezeType, 'InitFreeze') || strcmpi(freezeType, 'ManiFreeze')
            t_start_Early = freezEpisodes{frzi}.freezeTPhaseS(1) + 2;
        end
        t_end_Early = t_start_Early + 1;
        
        % Late Freeze
        t_end_Late = freezEpisodes{frzi}.freezeTPhaseS(2);
        t_start_Late = t_end_Late - 1;
        
        % middleFreeze
        t_start_Middle = t_end_Early;
        t_end_Middle = t_start_Late;


        % afterFreeze 200 ms 
        t_start_afterFreeze1s = t_end_Late;
        t_end_afterFreeze200ms = t_start_afterFreeze1s + 0.2;
        

        
        %%% segment freeze lfpdata into short leg: shortlfp

        % early Freeze segments
        t_str = t_start_Early;
        t_end = t_end_Early;
        idx_str = round(t_str * fs_lfp);
        idx_end = round(t_end * fs_lfp);
        seglfp  = lfpdata{tri}(:, idx_str: idx_end);
        
        shortlfp = [];
        len = size(seglfp, 2);
        shortSegn = floor(len / nseg);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nseg + 1;
            endi = shortsegi * nseg;
            shortlfp = cat(3, shortlfp, seglfp(:, stri : endi));
            clear stri endi
        end
        
        lfpsegs_freeze.(freezeType).earlyFreeze = cat(3, lfpsegs_freeze.(freezeType).earlyFreeze, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn
        
        
        % middle Freeze segments
        t_str = t_start_Middle;
        t_end = t_end_Middle;
        idx_str = round(t_str * fs_lfp);
        idx_end = round(t_end * fs_lfp);
        seglfp  = lfpdata{tri}(:, idx_str: idx_end);
        
        shortlfp = [];
        len = size(seglfp, 2);
        shortSegn = floor(len / nseg);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nseg + 1;
            endi = shortsegi * nseg;
            shortlfp = cat(3, shortlfp, seglfp(:, stri : endi));
            clear stri endi
        end
        
        lfpsegs_freeze.(freezeType).middleFreeze = cat(3, lfpsegs_freeze.(freezeType).middleFreeze, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn
        
        
        
        % late Freeze segments
        t_str = t_start_Late;
        t_end = t_end_Late;
        idx_str = round(t_str * fs_lfp);
        idx_end = round(t_end * fs_lfp);
        seglfp  = lfpdata{tri}(:, idx_str: idx_end);
        
        shortlfp = [];
        len = size(seglfp, 2);
        shortSegn = floor(len / nseg);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nseg + 1;
            endi = shortsegi * nseg;
            shortlfp = cat(3, shortlfp, seglfp(:, stri : endi));
            clear stri endi
        end
        
        lfpsegs_freeze.(freezeType).lateFreeze = cat(3, lfpsegs_freeze.(freezeType).lateFreeze, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn


        % after Freeze 200ms segments
        t_str = t_start_afterFreeze1s;
        t_end = t_end_afterFreeze200ms;
        idx_str = round(t_str * fs_lfp);
        idx_end = round(t_end * fs_lfp);
        seglfp  = lfpdata{tri}(:, idx_str: idx_end);
        
        shortlfp = [];
        len = size(seglfp, 2);
        shortSegn = floor(len / nseg);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nseg + 1;
            endi = shortsegi * nseg;
            shortlfp = cat(3, shortlfp, seglfp(:, stri : endi));
            clear stri endi
        end
        
        lfpsegs_freeze.(freezeType).afterFreeze200ms = cat(3, lfpsegs_freeze.(freezeType).afterFreeze200ms, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn

    end   
    
    clear nseg
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;
end


