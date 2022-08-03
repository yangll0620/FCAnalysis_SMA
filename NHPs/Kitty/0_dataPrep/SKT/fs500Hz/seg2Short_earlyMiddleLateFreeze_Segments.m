function [lfpsegs_freeze, fs, T_chnsarea, FreeTypes]= seg2Short_beforeEarlyMiddleLateFreeze_Segments(files, varargin)
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
if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpsegs_freeze = [];
    fs = [];
    T_chnsarea = [];
    FreeTypes = [];
    
    return;
end

% parse params
p = inputParser;
addParameter(p, 'tseg', 0.2, @isscalar);

parse(p,varargin{:});
tseg = p.Results.tseg;


FreeTypes = {'InitFreeze', 'ReachFreeze', 'ManipuFreeze'}; % combined {'freeze during React-Reach'}  and  {'freeze during Reach'} 

lfpsegs_freeze = struct();
lfpsegs_freeze.earlyFreeze = [];
lfpsegs_freeze.middleFreeze = [];
lfpsegs_freeze.lateFreeze = [];

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
        InitFreeze = false;
        ManiFreeze = true;
        if contains(freezeType, 'init Move') % 
            InitFreeze = true;
        end
        if contains(freezeType, 'Manipulation') %
            ManiFreeze = true;
        end
        
        
        % Early Freeze
        t_start_Early = freezEpisodes{frzi}.freezeTPhaseS(1);
        if InitFreeze || ManiFreeze
            t_start_Early = t_start_Early + 2;
        end
        t_end_Early = t_start_Early + 1;
        
        % Late Freeze
        t_end_Late = freezEpisodes{frzi}.freezeTPhaseS(2);
        t_start_Late = t_end_Late - 1;
        
        % middleFreeze
        t_end_Middle = freezEpisodes{frzi}.freezeTPhaseS(2) -1;
        t_start_Middle = freezEpisodes{frzi}.freezeTPhaseS(1) + 1;
        if InitFreeze || ManiFreeze
            t_start_Middle = t_start_Middle + 2;
        end
        
        
        %%% segment freeze lfpdata into short leg: shortlfp
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

    end   
    
    clear nseg
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;
end