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
files = dir(fullfile(inputfolder_FreezeEpisode, '*moderate*.mat'));


[lfpsegs_freeze, fs, T_chnsarea]= seg2Short_Freeze_Segments(files, 'tseg', seg_tseg);
savefilename = [savefile_prefix '-allFreeze.mat'];
savefile = fullfile(savefolder, savefilename);
save(savefile, 'lfpsegs_freeze', 'fs', 'T_chnsarea', 'seg_tseg');
clear lfpsegs_freeze fs T_chnsarea savefilename savefile


% lfpsegments of freeze phases using only reach freeze trials
[lfpsegs_reachFreeze, fs, T_chnsarea]= seg2Short_reachFreezeTrials_Segments(files, 'tseg', seg_tseg);
savefilename = [savefile_prefix '-reachFreezeTrials.mat'];
savefile = fullfile(savefolder, savefilename);
save(savefile, 'lfpsegs_reachFreeze', 'fs', 'T_chnsarea', 'seg_tseg');
clear lfpsegs_freeze fs T_chnsarea savefilename savefile


end


function [lfpsegs_freeze, fs, T_chnsarea]= seg2Short_reachFreezeTrials_Segments(files, varargin)
%
%  Just using reach Freeze trials to extract lfp segments
%
%   Input:
%       files: 
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
%             afterFreeze1s: [3×100×380 double]
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
lfpsegs_freeze.earlyFreeze = [];
lfpsegs_freeze.middleFreeze = [];
lfpsegs_freeze.lateFreeze = [];
lfpsegs_freeze.afterFreeze1s = [];

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

        freezeType = freezEpisodes{frzi}.freezeType;
        if contains(freezeType, 'init Move') || contains(freezeType, 'Manipulation')
            clear tri freezeType
            continue;
        end
        
        
        % Early Freeze
        t_start_Early = freezEpisodes{frzi}.freezeTPhaseS(1);
        t_end_Early = t_start_Early + 1;
        
        % Late Freeze
        t_end_Late = freezEpisodes{frzi}.freezeTPhaseS(2);
        t_start_Late = t_end_Late - 1;
        
        % middleFreeze
        t_start_Middle = t_end_Early;
        t_end_Middle = t_start_Late;


        % afterFreeze 1s 
        t_start_afterFreeze1s = t_end_Late;
        t_end_afterFreeze1s = t_start_afterFreeze1s + 1;
        

        
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
        
        lfpsegs_freeze.earlyFreeze = cat(3, lfpsegs_freeze.earlyFreeze, shortlfp);
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
        
        lfpsegs_freeze.middleFreeze = cat(3, lfpsegs_freeze.middleFreeze, shortlfp);
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
        
        lfpsegs_freeze.lateFreeze = cat(3, lfpsegs_freeze.lateFreeze, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn


        % after Freeze 1s segments
        t_str = t_start_afterFreeze1s;
        t_end = t_end_afterFreeze1s;
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
        
        lfpsegs_freeze.afterFreeze1s = cat(3, lfpsegs_freeze.afterFreeze1s, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn

    end   
    
    clear nseg
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;
end





function [lfpsegs_freeze, fs, T_chnsarea]= seg2Short_Freeze_Segments(files, varargin)
%
% 
%
%   Input:
%       files: 
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
%             afterFreeze1s: [3×100×380 double]
%            
%       fs, T_chnsarea
%

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
lfpsegs_freeze.earlyFreeze = [];
lfpsegs_freeze.middleFreeze = [];
lfpsegs_freeze.lateFreeze = [];
lfpsegs_freeze.afterFreeze1s = [];

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

        freezeType = freezEpisodes{frzi}.freezeType;
        
        
        % Early Freeze
        t_start_Early = freezEpisodes{frzi}.freezeTPhaseS(1);
        if contains(freezeType, 'init Move') || contains(freezeType, 'Manipulation')
            t_start_Early = t_start_Early + 2;
        end
        t_end_Early = t_start_Early + 1;
        
        % Late Freeze
        t_end_Late = freezEpisodes{frzi}.freezeTPhaseS(2);
        t_start_Late = t_end_Late - 1;
        
        % middleFreeze
        t_start_Middle = t_end_Early;
        t_end_Middle = t_start_Late;


        % afterFreeze 1s 
        t_start_afterFreeze1s = t_end_Late;
        t_end_afterFreeze1s = t_start_afterFreeze1s + 1;
        

        
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
        
        lfpsegs_freeze.earlyFreeze = cat(3, lfpsegs_freeze.earlyFreeze, shortlfp);
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
        
        lfpsegs_freeze.middleFreeze = cat(3, lfpsegs_freeze.middleFreeze, shortlfp);
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
        
        lfpsegs_freeze.lateFreeze = cat(3, lfpsegs_freeze.lateFreeze, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn


        % after Freeze 1s segments
        t_str = t_start_afterFreeze1s;
        t_end = t_end_afterFreeze1s;
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
        
        lfpsegs_freeze.afterFreeze1s = cat(3, lfpsegs_freeze.afterFreeze1s, shortlfp);
        clear t_str t_end idx_str idx_end seglfp shortlfp len shortSegn

    end   
    
    clear nseg
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;
end


