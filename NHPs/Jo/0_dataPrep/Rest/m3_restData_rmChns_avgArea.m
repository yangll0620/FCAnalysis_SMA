function m3_restData_rmChns_avgArea()
    %%  average lfp acros each area
    %
    %   1. replace STN into stn0-1, stn1-2 .. and GP into gp0-1....
    %   2. remove unwanted chns 
    %   3. average across areas
    %
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
    
    addpath(genpath(fullfile('.', 'subFuncs')));

    % the corresponding pipeline and the parent folder for this code
    [codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

    %% global variables
    % animal
    if ismac
        % Code to run on Mac platform
    elseif isunix
        % Code to run on Linux platform
        
        [fi, j] = regexp(codecorresfolder, ['NHPs', '/', '[A-Za-z]*']);
    elseif ispc
        % Code to run on Windows platform
        
        [fi, j] = regexp(codecorresfolder, ['NHPs', '\\', '[A-Za-z]*']);
    else
        disp('Platform not supported')
    end
    animal = codecorresfolder(fi + length('NHPs') + 1:j);

    %%  input setup

    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_Power');
    
    unwanted_DBS = unwanted_DBS_extract(animal);
    
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
        [data_segments, T_chnsarea_new] = rmChns_avgArea_1file(data_segments, T_chnsarea, unwanted_DBS);
        

        
        
        % save
        T_chnsarea = T_chnsarea_new;
        
        [i, j]= regexp(filename, [animal '_[a-zA-Z]*_' ]);
        savefilename = [animal '_' savefilename_addstr filename(j:end)];
        
        save(fullfile(savefolder, savefilename), 'data_segments',  'T_chnsarea', 'fs');
        
        
        clear filename file lfpdata T_chnsarea fs idx tempn savefilename
    end
end


function [avg_segments, T_chnsarea_new] = rmChns_avgArea_1file(data_segments, T_chnsarea, unWAreas)
    %%   remove the unwanted channels and average lfp across each area for data in file
    %
    % Args:
    %   data_segments: a struct storing data segment
    %   T_chnsarea: height = nchns
    %   unWAreas: unWanted areas e.g unWAreas = {'lCd', 'rMC', 'stn0-1', 'stn0-1', 'stn1-2', 'stn2-3', 'stn6-7', 'gp0-1'};
    %
    %  Returns:
    %   avg_segments: averaged lfp still struct (nsegs)
    %   T_chnsarea_new: new T_chnsarea_new (height = (nareas + nDBS))

    
    
    unWAreas_notDBS = unWAreas(~contains(unWAreas, 'stn') & ~contains(unWAreas, 'gp'));

    % extract mask_usedAreas (del unwanted areas, rm empty and DBS)
    mask_notDBSAreas = strcmp(T_chnsarea.electype, 'Gray Matter') | strcmp(T_chnsarea.electype, 'Utah Array');
    mask_unwanted = false(size(T_chnsarea.brainarea));
    for i = 1 : length(unWAreas_notDBS)
        uArea = unWAreas_notDBS{i};
        mask_unwanted = mask_unwanted | strcmp(T_chnsarea.brainarea, uArea);
        clear uArea
    end
    mask_emptyarea = cellfun(@(x) isempty(x), T_chnsarea.brainarea);
    mask_usedAreas = mask_notDBSAreas & ~mask_unwanted & ~mask_emptyarea;
    clear mask_notDBSAreas mask_unwanted mask_emptyarea unWAreas_DBS
    
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
    
    
    T_usedAreas = T_chnsarea(mask_usedAreas, :);
    T_DBS = T_chnsarea(mask_DBS, :);
    


    % extract the averaged lfp avglfp: (nareas + nDBS) * ntemp * nsegs
    avg_segments = struct(); 
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
        
        
        % assign
        avg_segments(segi).lfp = avglfp;
        
        clear lfp_1seg lfpdata_areas lfpdata_DBS avglfpdata_areas T_new
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
end


