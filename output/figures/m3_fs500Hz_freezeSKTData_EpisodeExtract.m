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



