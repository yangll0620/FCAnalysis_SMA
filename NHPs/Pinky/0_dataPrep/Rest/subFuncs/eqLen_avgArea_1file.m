function [eqLen_avglfp, T_chnsarea_new] = eqLen_avgArea_1file(data_segments, T_chnsarea, segt, fs)
    %%   segment into intervals with same length and average lfp across each area for data in file
    %
    % Args:
    %   data_segments: a struct storing data segment
    %   T_chnsarea: height = nchns
    %   segt: the same time length, unit second (e.g 2)
    %   fs: sample rate
    %
    %  Returns:
    %   eqLen_avglfp: averaged lfp with same length ((nareas + nDBS) * ntemp * nsegs)
    %   T_chnsarea_new: new T_chnsarea_new (height = (nareas + nDBS))


    % extract mask_usedAreas (del unwanted areas, rm empty and DBS)
    mask_notDBSAreas = strcmp(T_chnsarea.electype, 'Gray Matter') | strcmp(T_chnsarea.electype, 'Utah Array');
    mask_unwanted = strcmp(T_chnsarea.brainarea, 'lCd') | strcmp(T_chnsarea.brainarea, 'rMC');
    mask_emptyarea = cellfun(@(x) isempty(x), T_chnsarea.brainarea);
    mask_usedAreas = mask_notDBSAreas & ~mask_unwanted & ~mask_emptyarea;
    clear mask_notDBSAreas mask_unwanted mask_emptyarea
    
    
    mask_DBS = strcmp(T_chnsarea.electype, 'DBS');
    
    
    T_usedAreas = T_chnsarea(mask_usedAreas, :);
    T_DBS = T_chnsarea(mask_DBS, :);
    
    % the same length 
    ntemp_same = ceil(fs * segt);

    % extract the equal length averaged lfp eqLen_avglfp: (nareas + nDBS) * ntemp * nsegs
    eqLen_avglfp = []; 
    for segi = 1:length(data_segments)
        
        lfp_1seg = data_segments(segi).lfp; % lfp_1seg: ntemp * nchns 
        lfp_1seg = lfp_1seg'; % lfp_1seg: nchns * ntemp 
        
        
        lfpdata_areas = lfp_1seg(mask_usedAreas, :);
        lfpdata_DBS = lfp_1seg(mask_DBS, :);
        
        
        %%% avg lfp across each area, avglfpdata_areas: nareas * ntemp %%%%
        [avglfpdata_areas, T_new] = avglfp_acrossArea(lfpdata_areas, T_usedAreas); 
        
        if ~exist('T_usedAreas_new', 'var')
            T_usedAreas_new = T_new;
        else
            if ~(isequal(table2struct(T_new), table2struct(T_usedAreas_new)))
                disp(['T_new in' num2str(segi) ' not equal exist T_usedAreas_new.'])
                disp(T_new)
                disp(T_usedAreas_new)
            end
        end
        
        %%% cat notDBS areas and DBS %%%
        avglfp = cat(1, avglfpdata_areas, lfpdata_DBS);
        if segi == 1
            T_chnsarea_new = [T_usedAreas_new; T_DBS];
            T_chnsarea_new.chni = [1: height(T_chnsarea_new)]';
        end
        
        
        
        %%% seg into intervals with same time length   %%%%
        for segi_same = 1:floor(size(avglfpdata_areas, 2) / ntemp_same)
            idx_str = (segi_same - 1) * ntemp_same + 1;
            idx_end = segi_same * ntemp_same;

            eqLen_avglfp = cat(3, eqLen_avglfp, avglfp(:, idx_str:idx_end));

            clear idx_str idx_end
        end
             
        
        clear lfp_1seg lfpdata_areas lfpdata_DBS avglfpdata_areas T_new
    end
     
end