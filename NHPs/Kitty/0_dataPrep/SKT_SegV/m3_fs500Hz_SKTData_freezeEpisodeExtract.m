function m3_fs500Hz_SKTData_freezeEpisodeExtract()
% Objective:
%       extract freezing episode
%       data structure
%           episodes{}
% 1. manupulation time > 5s

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);



%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI_goodReach');


optFreezeTypes = {'freeze during init Move', 'freeze during React-Reach', 'freeze during Reach', 'freeze during Manipulation'};

pdcond = 'moderate';

speedThres_Move = 30;
tThes_Reaction = 5;
tThesFreeze_mani  = 5;


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code Start Here

% description stored inside and outside
discrip_out = ['optFreezeTypes: optional ' num2str(length(optFreezeTypes)) ' freezing Types; ' ...
    'tThesFreeze_mani: time threshold used for extracting freezing episode during manipulation (if larger than theshold, freezing episode)'];
discrip_inside = ['triali:trial number; ' ...
    'freezeType:current freezing type; ' ...
    'freezeTPhaseS:freezing time start and end time in second, idx = 1 is time 0;'];


files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
for fi = 1 : length(files)
    filename = files(fi).name;
    
    load(fullfile(inputfolder, filename), 'lfpdata', 'T_chnsarea', 'T_idxevent_lfp', 'fs_lfp',  ...
        'Wpos_smooth_trial', 'Wrist_smooth_trial', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma', ...
        'selectedTrials');
    
    
    
    %%% init freezStruct
    freezStruct = struct();
    freezStruct.optFreezeTypes = optFreezeTypes;
    freezStruct.tThesFreeze_mani = tThesFreeze_mani;
    freezStruct.discription = discrip_out;
    
    
    freezEpisodes = {};
    
    
    
    %%% --- extract freezing episodes during reaction phase
    times_reaction = (T_idxevent_ma.ReachTimeix - T_idxevent_ma.TargetTimeix) / fs_ma;
    idxs_freeze_reaction = find(times_reaction >= tThes_Reaction);
    for idi = 1 : length(idxs_freeze_reaction)
        tri = idxs_freeze_reaction(idi);
        
        
        %%% --- find the init move freezing and reach freezing seperately -- %%%
        
        % extract speed vector during reaction phase
        speeds_inReaction = smoothWspeed_trial{tri}(T_idxevent_ma.TargetTimeix(tri):T_idxevent_ma.ReachTimeix(tri), 1);
        idxs_move = find(speeds_inReaction >= speedThres_Move); % idxs for moving
        if isempty(idxs_move)
            idxs_StrEnd_freeze = [1 length(speeds_inReaction)];
        else
            % extract idxs_StrEnd_moves
            diffIdxs = [0; diff(idxs_move)];
            idxs_str = find(diffIdxs ~=1); % move start index for each move phase
            if length(idxs_str) > 1
                idxs_end = [idxs_str(2:end) - 1; length(idxs_move)];
            else
                idxs_end = length(idxs_move);
            end
            idxs_StrEnd_moves = [idxs_move(idxs_str), idxs_move(idxs_end)];
            
            % extract idxs_StrEnd_freeze
            idxs_StrEnd_freeze = [];
            for idi = 1 : size(idxs_StrEnd_moves, 1)
                if idi == 1 &&  idxs_StrEnd_moves(idi, 1) ~=1 % first one
                    idxs_StrEnd_freeze = [idxs_StrEnd_freeze;1 idxs_StrEnd_moves(idi, 1)-1];
                end
                
                if idi > 1
                    idxs_StrEnd_freeze = [idxs_StrEnd_freeze;idxs_StrEnd_moves(idi-1, 2)+1 idxs_StrEnd_moves(idi, 1)-1];
                end
                
                if idi == size(idxs_StrEnd_moves, 1) && idxs_StrEnd_moves(idi, 2) ~= length(speeds_inReaction) % last one
                    idxs_StrEnd_freeze = [idxs_StrEnd_freeze;idxs_StrEnd_moves(idi, 2)+1 length(speeds_inReaction)];
                end
            end
            clear diffIdxs idxs_str idxs_end idxs_StrEnd_moves
        end
        clear speeds_inReaction idxs_move
        
        % freeze during init move
        s = struct();
        freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'init Move', 'IgnoreCase',true), optFreezeTypes)};
        s.freezeType = freezeType;
        s.triali = tri;
        s.filename = filename;
        s.folder = inputfolder;
        s.freezeTPhaseS = idxs_StrEnd_freeze(1, :) / fs_ma;
        s.discription = discrip_inside;
        freezEpisodes{end + 1, 1} = s;
        clear s freezeType 

        
        % freeze during reach
        if size(idxs_StrEnd_freeze, 1) > 1
            s = struct();
            freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'React-Reach', 'IgnoreCase',true), optFreezeTypes)};
            s.freezeType = freezeType;
            s.triali = tri;
            s.filename = filename;
            s.folder = inputfolder;
            s.freezeTPhaseS = idxs_StrEnd_freeze(2:end, :) / fs_ma;
            s.discription = discrip_inside;
            freezEpisodes{end + 1, 1} = s;
            clear s freezeType 
        end
        
        clear idxs_StrEnd_freeze tri
    end
    clear idxs_freeze_reaction times_reaction
    
    
    %%% --- extract freezing episodes during manupation phase
    freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'manipulation', 'IgnoreCase',true), optFreezeTypes)};
    times_mani = (T_idxevent_ma.ReturnTimeix - T_idxevent_ma.TouchTimeix) / fs_ma;
    idxs_freeze_mani = find(times_mani >= tThesFreeze_mani);
    for idi = 1 : length(idxs_freeze_mani)
        tri = idxs_freeze_mani(idi);
        
        s = struct();
        
        s.freezeType = freezeType;
        s.triali = tri;
        s.filename = filename;
        s.folder = inputfolder;
        s.freezeTPhaseS = [T_idxevent_ma.TouchTimeix(tri) T_idxevent_ma.ReturnTimeix(tri)] / fs_ma;
        s.discription = discrip_inside;
        
        freezEpisodes{end + 1, 1} = s;
        
        clear s tri
    end
    clear freezeType times_mani idxs_freeze_mani
    
    %%% assign to freezStruct.freezEpisodes
    freezStruct.freezEpisodes = freezEpisodes;
    
    tmp = regexp(filename, '\d{8}', 'match');
    dateofexp_str = tmp{1};
    tmp = regexp(filename, 'bktdt\d{1}', 'match');
    bktdt_str = tmp{1};
    clear tmp
    
    save(fullfile(savefolder, [animal '_freezeEpisodes_' pdcond '_' dateofexp_str '_' bktdt_str '.mat']), 'lfpdata', 'T_chnsarea', 'T_idxevent_lfp', 'fs_lfp',  ...
        'Wpos_smooth_trial', 'Wrist_smooth_trial', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma', ...
        'selectedTrials', 'freezStruct');
    
    clear('lfpdata', 'T_chnsarea', 'T_idxevent_lfp', 'fs_lfp',  ...
        'Wpos_smooth_trial', 'Wrist_smooth_trial', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma', ...
        'selectedTrials');
    
    clear filename freezStruct dateofexp_str bktdt_str
    
end





