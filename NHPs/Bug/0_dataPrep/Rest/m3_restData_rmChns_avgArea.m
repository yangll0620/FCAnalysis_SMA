function m3_restData_rmChns_avgArea()
    %%  average lfp acros each area
    %
    %   1. combine Vlo and VPLo
    %   2. replace STN into stn0-1, stn1-2 .. and GP into gp0-1....
    %   3. average across areas (M1 and PMC in selected Layer)
    %
    % Input:
    %   m2_restData_selectSeg_Power
    %
    % Outputs:
    %   lfpdata: averaged lfp with same length ((nareas + nDBS) * ntemp * nsegs)
    %   T_chnsarea: new T_chnsarea_new (height = (nareas + nDBS))
    %   fs: sample rate 


    %% folders generate
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');

    % find the codefolder
    idx = strfind(codefilepath, 'code');
    codefolder = codefilepath(1:idx + length('code') - 1);
    clear idx

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));


    % the corresponding pipeline and the parent folder for this code
    [codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

    %% global variables
    % animal
    [fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
    animal = codecorresfolder(fi + length('NHPs/'):j);

    segt = 2;

    %%  input setup

    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_Power');
    
    unWAreas = {};

    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = 'selAreas_avgArea';

    %% starting: narrow filter the lfp data of all the files
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);
    
    
    for filei = 1:nfiles
        
        filename = files(filei).name;
        file = fullfile(files(filei).folder, files(filei).name);
           
        
        load(file, 'data_segments', 'fs', 'T_chnsarea')
        
        if isempty(data_segments)

            continue;
        end
        
        
        % replace VLo and VPLo with VLoVPLo
        mask_VLoVPLo = strcmp(T_chnsarea.brainarea, 'VLo') | strcmp(T_chnsarea.brainarea, 'VPLo');
        T_chnsarea{mask_VLoVPLo, 'brainarea'} = {'VLo/VPLo'};

        
        
        % replace STN into stn0-1, stn1-2 ... 
        row_STN = find(strcmp(T_chnsarea.brainarea, 'STN'));
        for i = 1: length(row_STN)
            recChn = T_chnsarea.recordingchn(row_STN(i));
            T_chnsarea.brainarea(row_STN(i)) = {['stn' num2str(recChn-1) '-' num2str(recChn)]};
            clear recChn
        end
        
        % replace GP into gp0-1, gp1-2 ... 
        row_STN = find(strcmp(T_chnsarea.brainarea, 'GP'));
        for i = 1: length(row_STN)
            recChn = T_chnsarea.recordingchn(row_STN(i));
            T_chnsarea.brainarea(row_STN(i)) = {['gp' num2str(recChn-1) '-' num2str(recChn)]};
            clear recChn
        end
        
        
        disp(filename);
        
        % rm unwanted areas and average across areas
        [lfpdata, T_chnsarea_new] = rmChns_avgArea_1file(data_segments, T_chnsarea, unWAreas);
        if isempty(lfpdata)
            continue;
        end
        

        
        
        % save
        T_chnsarea = T_chnsarea_new;
        
        idx = strfind(filename, [animal '_']);
        tmpn = length([animal '_']);
        savefilename = [filename(idx:idx + tmpn - 1) savefilename_addstr ...
            upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];
        
        save(fullfile(savefolder, savefilename), 'lfpdata',  'T_chnsarea', 'fs');
        
        
        clear filename file lfpdata T_chnsarea fs idx tempn savefilename
    end
end


function [avg_segments, T_chnsarea_new] = rmChns_avgArea_1file(data_segments, T_chnsarea, unWAreas)
    %%   remove the unwanted channels and average lfp across each area for data in file
    %
    % Args:
    %   data_segments: a struct storing data segment
    %   T_chnsarea: height = nchns
    %   
    %
    %  Returns:
    %   avg_segments: averaged lfp still struct (nsegs), selected Layer useds for M1 and PMC
    %   T_chnsarea_new: new T_chnsarea_new (height = (nareas + nDBS))
    
    
    depth_M1Layer5 = [10 14];
    depth_PMCLayer5 = [10 14];
    
    depthUsed_M1 = depth_M1Layer5;
    depthUsed_PMC = depth_PMCLayer5;
    
    
    unWAreas_notDBS = unWAreas(~contains(unWAreas, 'stn') & ~contains(unWAreas, 'gp'));

    % extract mask_usedAreas (del unwanted areas, rm empty and DBS)
    mask_notDBSAreas = ~strcmp(T_chnsarea.electype, 'DBS');
    mask_unwanted = false(size(T_chnsarea.brainarea));
    for i = 1 : length(unWAreas_notDBS)
        uArea = unWAreas_notDBS{i};
        mask_unwanted = mask_unwanted | strcmp(T_chnsarea.brainarea, uArea);
        clear uArea
    end
    mask_emptyarea = cellfun(@(x) isempty(x), T_chnsarea.brainarea);
    mask_M1PMC = strcmp(T_chnsarea.brainarea, 'M1') | strcmp(T_chnsarea.brainarea, 'PMC');
    mask_usedAreas = mask_notDBSAreas & ~mask_unwanted & ~mask_emptyarea;
    mask_usedAreas_notM1PMC = mask_usedAreas & ~mask_M1PMC;
    mask_usedAreas_M1PMC = mask_usedAreas & mask_M1PMC;
    clear mask_notDBSAreas mask_unwanted mask_emptyarea unWAreas_DBS mask_usedAreas
    
    % DBS (del unwanted stn contact)
    unWAreas_DBS = unWAreas(contains(unWAreas, 'stn') | contains(unWAreas, 'gp'));
    mask_unwanted = false(size(T_chnsarea.brainarea));
    for i = 1 : length(unWAreas_DBS)
        uArea = unWAreas_DBS{i};
        mask_unwanted = mask_unwanted | strcmp(T_chnsarea.brainarea, uArea);
        clear uArea
    end
    mask_DBS = strcmp(T_chnsarea.electype, 'DBS') & ~mask_unwanted;
    clear unWAreas_DBS mask_unwanted
    
    
    T_usedAreas_notM1PMC = T_chnsarea(mask_usedAreas_notM1PMC, :);
    T_usedAreas_M1PMC = T_chnsarea(mask_usedAreas_M1PMC, :);
    T_DBS = T_chnsarea(mask_DBS, :);
    


    % extract the averaged lfp avglfp seg
    avg_segments = struct(); 
    for segi = 1:length(data_segments)
        
        lfp_1seg = data_segments(segi).lfp; % lfp_1seg: ntemp * nchns 
        lfp_1seg = lfp_1seg'; % lfp_1seg: nchns * ntemp 
        
        
        lfpdata_areas_notM1PMC = lfp_1seg(mask_usedAreas_notM1PMC, :);
        lfpdata_areas_M1PMC = lfp_1seg(mask_usedAreas_M1PMC, :);
        lfpdata_DBS = lfp_1seg(mask_DBS, :);
        
        
        
        %%% avg lfp across each area (M1 and PMC), avglfpdata_areas: nareas * ntemp %%%%
        [avglfpdata_areas_M1PMC, T_new_M1PMC] = avglfp_acrossM1PMC(lfpdata_areas_M1PMC, T_usedAreas_M1PMC, depthUsed_M1, depthUsed_PMC);
        if isempty(avglfpdata_areas_M1PMC)
            disp('avglfpdata_areas_M1PMC is empty, return!')
            avg_segments = [];
            T_chnsarea_new = [];
            return ;
        end
        if segi == 1
            T_usedAreasM1PMC_new = T_new_M1PMC;
        else
            if ~(isequal(table2struct(T_new_M1PMC), table2struct(T_usedAreasM1PMC_new)))
                disp(['T_new in' num2str(segi) ' not equal exist T_usedAreas_new.'])
                disp(T_new_M1PMC)
                disp(T_usedAreasM1PMC_new)
            end
        end
         
        %%% avg lfp across each area, avglfpdata_areas: nareas * ntemp %%%%
        [avglfpdata_areas_notM1PMC, T_new_notM1PMC] = avglfp_acrossArea(lfpdata_areas_notM1PMC, T_usedAreas_notM1PMC); 
        
        if segi == 1
            T_usedAreasNotM1PMC_new = T_new_notM1PMC;
        else
            if ~(isequal(table2struct(T_new_notM1PMC), table2struct(T_usedAreasNotM1PMC_new)))
                disp(['T_new in' num2str(segi) ' not equal exist T_usedAreas_new.'])
                disp(T_new_notM1PMC)
                disp(T_usedAreasNotM1PMC_new)
            end
        end
        
            
        
        %%% cat notDBS M1PMC, notM1PMC areas and DBS %%%
        avglfp = cat(1, avglfpdata_areas_M1PMC, avglfpdata_areas_notM1PMC, lfpdata_DBS);
        if segi == 1
            T_chnsarea_new = [T_usedAreasM1PMC_new; T_usedAreasNotM1PMC_new; T_DBS];
            T_chnsarea_new.chni = [1: height(T_chnsarea_new)]';
        end
        
        
        % assign
        avg_segments(segi).lfp = avglfp;
        
        clear lfp_1seg lfpdata_areas* lfpdata_DBS avglfpdata_areas T_new*
    end
    

     
end


function [avglfp, T_chnsarea_new] = avglfp_acrossArea(lfpdata, T_chnsarea)
    % average lfp across each area 
    % 
    % args:
    %       lfpdata : lfp data in area (nchns * ntemp)
    %       T_chnsarea: chns-area table 
    % return:
    %       avglfp: averaged lfp data (nareas * ntemp)
    %       T_chnsarea_new: the new chns-area table


    uniqAreas = unique(T_chnsarea.brainarea);
    avglfp = [];
    T_chnsarea_new = T_chnsarea([], :);
    for ai = 1 : length(uniqAreas)
        barea = uniqAreas{ai};


        mask_area = strcmp(T_chnsarea.brainarea, barea);
        lfp_area = lfpdata(mask_area, :);
        avglfp_area = mean(lfp_area, 1);


        avglfp = cat(1, avglfp, avglfp_area);

        tbl = T_chnsarea(find(mask_area,1),:);
        tbl.recordingchn = nan;
        T_chnsarea_new = [T_chnsarea_new; tbl];

        clear barea mask_area avglfp_area tbl
    end
    T_chnsarea_new.chni = [1: height(T_chnsarea_new)]';
end % avglfp_acrossArea

function [avglfp, T_chnsareaM1PMC_new] = avglfp_acrossM1PMC(lfpdata, T_chnsareaM1PMC, depthUsed_M1, depthUsed_PMC)
    % average lfp across M1 or PMC based on depthUsed 
    % 
    % args:
    %       lfpdata : lfp data in area (nchns * ntemp)
    %       T_chnsarea: chns-area table 
    % return:
    %       avglfp: averaged lfp data (nareas * ntemp)
    %       T_chnsarea_new: the new chns-area table


    uniqAreas = unique(T_chnsareaM1PMC.brainarea);
    avglfp = [];
    T_chnsareaM1PMC_new = T_chnsareaM1PMC([], :);
    for ai = 1 : length(uniqAreas)
        barea = uniqAreas{ai};
        
        if ~strcmp(barea, 'M1') && ~strcmp(barea, 'PMC')
            avglfp = []; T_chnsareaM1PMC_new = [];
            disp([barea ' exist, only M1 and PMC are allowed!']);
            return;
        end
                
        if strcmp(barea, 'M1')
            depthUsed = depthUsed_M1 ;
        end
        if strcmp(barea, 'PMC')
            depthUsed = depthUsed_PMC ;
        end
        
        mask_area = strcmp(T_chnsareaM1PMC.brainarea, barea) & T_chnsareaM1PMC.depth >= depthUsed(1) & T_chnsareaM1PMC.depth <= depthUsed(2);
        if ~any(mask_area) % if not channels in depthUsed
            avglfp = []; T_chnsareaM1PMC_new = [];
            disp([ 'no channels in Useddept for  ' barea]);
            return;
        end

        lfp_area = lfpdata(mask_area, :);
        avglfp_area = mean(lfp_area, 1);


        avglfp = cat(1, avglfp, avglfp_area);

        tbl = T_chnsareaM1PMC(find(mask_area,1),:);
        tbl.recordingchn = nan;
        T_chnsareaM1PMC_new = [T_chnsareaM1PMC_new; tbl];

        clear barea mask_area avglfp_area tbl
    end
    T_chnsareaM1PMC_new.chni = [1: height(T_chnsareaM1PMC_new)]';
end






