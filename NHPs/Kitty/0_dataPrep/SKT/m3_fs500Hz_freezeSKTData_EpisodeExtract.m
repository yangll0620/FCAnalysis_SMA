function m3_fs500Hz_freezeSKTData_EpisodeExtract()
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


pdcond = 'moderate';

speedThres_Move = 30;

tThesFreeze_init = 5;
tThesFreeze_reach  = 3;
tThesFreeze_mani  = 5;



%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
saveFreezTrialsfolder = fullfile(savefolder, ['freezedTrials-tThesFreezeReach' num2str(tThesFreeze_reach) 's']);
copyfile2folder(codefilepath, savecodefolder);
if ~exist(saveFreezTrialsfolder, 'dir')
    mkdir(saveFreezTrialsfolder);
end
savefilename_prefix = [animal '_freezeEpisodes_' pdcond '-tThesFreezeReach' num2str(tThesFreeze_reach) 's'];

%% Code Start Here

optFreezeTypes = optFreezeTypes_extract('codesavefolder', savecodefolder);

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
    freezStruct.speedThres_Move = speedThres_Move;
    freezStruct.tThesFreeze_init = tThesFreeze_init;
    freezStruct.tThesFreeze_reach = tThesFreeze_reach;
    freezStruct.tThesFreeze_mani = tThesFreeze_mani;
    freezStruct.discription = discrip_out;
    
    
    %%% extract freezEpisodes %%%
    freezEpisodes = {};
    for tri = 1: length(smoothWspeed_trial)    
        
        %%% --- extract freezing episodes during reaction phase: initFreeze and reachFreeze
        idx_strInTrial = T_idxevent_ma.TargetTimeix(tri);
        idx_endInTrial = T_idxevent_ma.ReachTimeix(tri);
        speeds_inReaction = smoothWspeed_trial{tri}(idx_strInTrial:idx_endInTrial, 1);
        idxs_nomove = find(speeds_inReaction < speedThres_Move); % idxs for not moving
        
        if ~isempty(idxs_nomove)
            % seprate into nomove intervals
            diffIdxs = [0; diff(idxs_nomove)];
            idxs_str = find(diffIdxs ~=1);
            if length(idxs_str) > 1
                idxs_end = [idxs_str(2:end) - 1; length(idxs_nomove)];
            else
                idxs_end = length(idxs_nomove);
            end
            idxs_StrEnd_moves = [idxs_nomove(idxs_str), idxs_nomove(idxs_end)];
            
            % append freeze episode, separate into initFreeze and reachFreeze
            for idxi = 1 : size(idxs_StrEnd_moves, 1)
                t_freeze = (idxs_StrEnd_moves(idxi, 2) - idxs_StrEnd_moves(idxi, 1)) / fs_ma;
                
                % separate into initFreeze and reachFreeze
                if idxi == 1 && idxs_StrEnd_moves(idxi, 1) == 1
                    if t_freeze >= tThesFreeze_init % init Freeze case
                        freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'init Move', 'IgnoreCase',true), optFreezeTypes)};
                        
                        s = struct();
                        s.freezeType = freezeType;
                        s.triali = tri;
                        s.filename = filename;
                        s.folder = inputfolder;
                        s.freezeTPhaseS = (idxs_StrEnd_moves(idxi, :) + idx_strInTrial -1)/ fs_ma;
                        s.discription = discrip_inside;
                        
                        freezEpisodes{end + 1, 1} = s;
                        
                        clear freezeType s
                    end
                    
                else
                    if t_freeze >= tThesFreeze_reach % reach freeze in reaction phase case
                        freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'freeze during React-Reach', 'IgnoreCase',true), optFreezeTypes)};
                        
                        s = struct();
                        s.freezeType = freezeType;
                        s.triali = tri;
                        s.filename = filename;
                        s.folder = inputfolder;
                        s.freezeTPhaseS = (idxs_StrEnd_moves(idxi, :) + idx_strInTrial -1)/ fs_ma;
                        s.discription = discrip_inside;
                        
                        freezEpisodes{end + 1, 1} = s;
                        
                        clear freezeType s
                    end
                    clear t_freeze
                end
                clear t_freeze
            end
            clear diffIdxs idxs_str idxs_end idxs_StrEnd_moves idxi
        end
        clear idx_strInTrial idx_endInTrial speeds_inReaction idxs_nomove
        
        %%% --- extract freezing episodes during Reach phase
        idx_strInTrial = T_idxevent_ma.ReachTimeix(tri);
        idx_endInTrial = T_idxevent_ma.TouchTimeix(tri);
        speeds_inReach = smoothWspeed_trial{tri}(idx_strInTrial:idx_endInTrial, 1);
        idxs_nomove = find(speeds_inReach < speedThres_Move); % idxs for not moving
        
        if ~isempty(idxs_nomove)
            % seprate into nomove intervals
            diffIdxs = [0; diff(idxs_nomove)];
            idxs_str = find(diffIdxs ~=1);
            if length(idxs_str) > 1
                idxs_end = [idxs_str(2:end) - 1; length(idxs_nomove)];
            else
                idxs_end = length(idxs_nomove);
            end
            idxs_StrEnd_moves = [idxs_nomove(idxs_str), idxs_nomove(idxs_end)];
            
            % append freeze during Reach episode
            for idxi = 1 : size(idxs_StrEnd_moves, 1)
                t_freeze = (idxs_StrEnd_moves(idxi, 2) - idxs_StrEnd_moves(idxi, 1)) / fs_ma;
                if t_freeze >= tThesFreeze_reach
                    freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'freeze during Reach', 'IgnoreCase',true), optFreezeTypes)};
                    
                    s = struct();
                    s.freezeType = freezeType;
                    s.triali = tri;
                    s.filename = filename;
                    s.folder = inputfolder;
                    s.freezeTPhaseS = (idxs_StrEnd_moves(idxi, :) + idx_strInTrial -1)/ fs_ma;
                    s.discription = discrip_inside;
                    
                    freezEpisodes{end + 1, 1} = s;
                    
                    clear freezeType s
                end
                clear t_freeze
            end
            clear diffIdxs idxs_str idxs_end idxs_StrEnd_moves idxi
        end
        clear idx_strInTrial idx_endInTrial speeds_inReach idxs_nomove
        
        
        %%% --- extract freezing episodes during manupation phase
        t_mani = (T_idxevent_ma.ReturnTimeix(tri) - T_idxevent_ma.TouchTimeix(tri)) / fs_ma;
        if t_mani >= tThesFreeze_mani % append manipulation freeze episode
            freezeType = optFreezeTypes{cellfun(@(x) contains(x, 'manipulation', 'IgnoreCase',true), optFreezeTypes)};
            
            s = struct();
            s.freezeType = freezeType;
            s.triali = tri;
            s.filename = filename;
            s.folder = inputfolder;
            s.freezeTPhaseS = [T_idxevent_ma.TouchTimeix(tri) T_idxevent_ma.ReturnTimeix(tri)]/ fs_ma;
            s.discription = discrip_inside;
            
            freezEpisodes{end + 1, 1} = s;
            
            clear s freezeType
        end
        clear t_mani
    end
    
    %%% assign to freezStruct.freezEpisodes
    freezStruct.freezEpisodes = freezEpisodes;
    clear freezEpisodes
    
    
    %%% save part %%%
    tmp = regexp(filename, '\d{8}', 'match');
    datebktdt_str = tmp{1};
    tmp = regexp(filename, 'bktdt\d{1}', 'match');
    bktdt_str = tmp{1};
    clear tmp
    
    save(fullfile(savefolder, [savefilename_prefix '_' datebktdt_str '_' bktdt_str '.mat']), 'lfpdata', 'T_chnsarea', 'T_idxevent_lfp', 'fs_lfp',  ...
        'Wpos_smooth_trial', 'Wrist_smooth_trial', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma', ...
        'selectedTrials', 'freezStruct');
    
    
    %%% final clear
    clear('lfpdata', 'T_chnsarea', 'T_idxevent_lfp', 'fs_lfp',  ...
        'Wpos_smooth_trial', 'Wrist_smooth_trial', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma', ...
        'selectedTrials');
    clear filename freezStruct dateofexp_str bktdt_str
    
end


%%% -- Plot and check some of the extracted Freezed trials --
eveNames = {'cueonset', 'reachonset', 'touch', 'returnonset', 'mouth'};
initfrezPhaNames = {'cueonset', 'initFreezeEnd'};
reachfrezPhaNames = {'reachFreezStart', 'reachFreezEnd'};
manifrezPhaNames = {'touch', 'maniFreezEnd'};
freezCols = {'k', 'b', 'r'};
frezfiles = dir(fullfile(savefolder, [savefilename_prefix '*.mat']));
for fi = 1 : length(frezfiles)
    filename = frezfiles(fi).name;
    load(fullfile(savefolder, filename), 'freezStruct', ...
        'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma', 'selectedTrials');
    
    datebktdt_str = regexp(filename, '\d{8}_bktdt\d{1}', 'match');
    datebktdt_str = strrep(datebktdt_str{1}, '_', '-');
    
    freezEpisodes = freezStruct.freezEpisodes;
    tri_pre = 0;
    
    for frzi = 1: length(freezEpisodes)
        tri = freezEpisodes{frzi}.triali;
        if ~selectedTrials(tri)
            clear tri
            continue;
        end
        
        if tri ~= tri_pre % new trial
            if tri_pre ~= 0 % save and close previous trial
                
                annotation(gcf,'textbox',[0.01 0.8 0.1 0.2], 'LineStyle','none', ...
                    'String',{['tThesFreeze-init = ' num2str(freezStruct.tThesFreeze_init) 's'], ...
                              ['tThesFreeze-reach = ' num2str(freezStruct.tThesFreeze_reach) 's'], ...
                              ['tThesFreeze-mani = ' num2str(freezStruct.tThesFreeze_mani) 's']});
                
                legend(hlegshows, 'Orientation', 'horizontal', 'Location', 'north')
                saveas(gcf, fullfile(saveFreezTrialsfolder, [animal '-' datebktdt_str '-trial' num2str(tri_pre) '.tif']));
                clear ax
                close gcf
                clear hlegshows reachLegExist
                clear frzi_InTrial  
            end
            
            tri_pre = tri;
            hlegshows = []; % store the objects whose legend to show
            reachLegExist = false;
            
            
            figure('Position', [50 150 1800 600]);
            ax = axes(gcf);
            
            % plot smoothWspeed_trial
            ma = smoothWspeed_trial{tri};
            hp = plot(ax, [1: length(ma)]/fs_ma, ma, 'DisplayName','speed'); hold on
            hlegshows = [hlegshows hp];
            hp = plot(xlim, [speedThres_Move speedThres_Move], 'b-.', 'DisplayName','speedThres');
            hlegshows = [hlegshows hp];
            clear ma hp
            
            
            % plot event line
            xtks = xticks(ax);
            xtklabs = xticklabels(ax);
            tevents_ma = T_idxevent_ma{tri_pre, :} / fs_ma;
            for tei = 1 : length(tevents_ma)
                tevent = tevents_ma(tei);
                plot([tevent tevent], ylim, '--');
                xtks = [xtks tevent];
                xtklabs = [xtklabs; eveNames{tei}];
                clear tevent
            end
            [xtks, idxs]= sort(xtks);
            xtklabs = xtklabs(idxs);
            set(ax,'XTick', xtks, 'XTickLabel', xtklabs, 'XTickLabelRotation', 45);
            clear xtks xtklabs tevents_ma tei
            
            % title
            title([animal '-' datebktdt_str ', tri = ' num2str(tri)])
            
            frzi_InTrial = 1;
        else
            frzi_InTrial = frzi_InTrial + 1;
        end
        
        
        %%% plot start and end freeze lines
        ys = ylim;
        freezeType = freezEpisodes{frzi}.freezeType;
        idx = find(cellfun(@(x) contains(freezeType, x), {'init', 'Reach', 'Manipulation'}));
        Tphase = freezEpisodes{frzi}.freezeTPhaseS;
        switch idx
            case 1
                frezPhaNames = initfrezPhaNames;
                ys = [ys(1) ys(2)*2/3];
                
                % plot lines
                hp = plot([Tphase(1) Tphase(1)], ys, [freezCols{idx} '-'], 'DisplayName', 'initFreeze');
                plot([Tphase(2) Tphase(2)], ys, [freezCols{idx} '-']);
                hlegshows = [hlegshows hp];
                
                clear hp
            case 2
                frezPhaNames = reachfrezPhaNames;
                ys = [ys(1) ys(2)*1/2];
                
                % plot lines
                hp = plot([Tphase(1) Tphase(1)], ys, [freezCols{idx} '-'], 'DisplayName', 'reachFreeze');
                plot([Tphase(2) Tphase(2)], ys, [freezCols{idx} '-']);
                if ~reachLegExist
                    hlegshows = [hlegshows hp];
                    reachLegExist = true;
                end
            case 3
                frezPhaNames = manifrezPhaNames;
                ys = [ys(1) ys(2)*1/3];
                
                % plot lines
                hp = plot([Tphase(1) Tphase(1)], ys, [freezCols{idx} '-'], 'DisplayName', 'maniFreeze');
                plot([Tphase(2) Tphase(2)], ys, [freezCols{idx} '-']);
                hlegshows = [hlegshows hp];
                
                clear hp
            otherwise
                disp(['freeze type is not right: ' freezeType])
        end
         
        % plot texts
        str = ['fStart' num2str(frzi_InTrial)];
        ht = text(ax, Tphase(1), ys(2), str);
        set(ht, 'Units', 'points')
        pos_points = get(ht, 'Position');
        set(ht, 'Position', [pos_points(1)-length(str)/2*5  pos_points(2)+10]);
        str = ['fEnd' num2str(frzi_InTrial)];
        ht = text(ax, Tphase(2), ys(2) + 2, str);
        set(ht, 'Units', 'points')
        pos_points = get(ht, 'Position');
        set(ht, 'Position', [pos_points(1)-length(str)/2*5  pos_points(2)+10]);
        clear ht str
        
        
        if frzi == length(freezEpisodes) % last frzi
            annotation(gcf,'textbox',[0.01 0.8 0.1 0.2], 'LineStyle','none', ...
                    'String',{['tThesFreeze-init = ' num2str(freezStruct.tThesFreeze_init) 's'], ...
                              ['tThesFreeze-reach = ' num2str(freezStruct.tThesFreeze_reach) 's'], ...
                              ['tThesFreeze-mani = ' num2str(freezStruct.tThesFreeze_mani) 's']});
            legend(hlegshows, 'Orientation', 'horizontal','Location', 'north')
            saveas(gcf, fullfile(saveFreezTrialsfolder, [animal '-' datebktdt_str '-trial' num2str(tri) '.tif']));
            clear ax
            close gcf
            clear hlegshows reachLegExist
            clear frzi_InTrial
        end
        
        clear frezPhaNames idx freezeType
    end
    
    clear filename freezStruct 
    clear('freezStruct', 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma');
    clear datebktdt_str 
    clear freezEpisodes tri_pre
end



