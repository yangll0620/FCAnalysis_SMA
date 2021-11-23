function m3_restData_rmChns_avgArea()
    %% extract m3_restData_rmChns_avgArea data directly using Ying's Kitty Rest data
    %
    %   1. combine avg_lfpm1 and bipolar lfpstn and lfpgp (data in Ying's file are already averaged lfpm1 and bipolar DBS data)
    %
    %   2. Down sample into 500Hs from 3000Hz
    %
    %   2. remove unwanted chns 
    %
    %   3. save based on individual dateBlock pair
    %
    %   Inputs:
    %       root2\Ying Yu\SKB_EventRelatedBeta\datasets_Allanimals\Resting\data_segment_Kitty_Rest_NonMovwithMAbyCleandataNo60HzFilt.mat
    %
    %   Outputs:
    %       lfpdata: averaged lfp with same length ((nareas + nDBS) * ntemp * nsegs)
    %       T_chnsarea: (height = (nareas + nDBS))
    %       fs: sample rate 
    
    
    
    
    %% folder generate
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');

    % find the codefolder
    idx = strfind(codefilepath, 'code');
    codefolder = codefilepath(1:idx + length('code') - 1);
    clear idx

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));
    addpath(genpath(fullfile(codefolder, 'NHPs')));

    % the corresponding pipeline folder for this code
    [codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

    % datafolder
    [datafolder, ~, ~, ~] = exp_subfolders();

    %% global parameter
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


    %% input setup

    % input folder
    input_folder = fullfile(datafolder, animal, 'root2', 'Ying Yu', 'SKB_EventRelatedBeta', 'datasets_Allanimals', 'Resting');
    input_file = 'data_segment_Kitty_Rest_NonMovwithMAbyCleandataNo60HzFilt.mat';
    

    unwanted_DBS = unwanted_DBS_extract(animal);
        
    fs_infile = 3000;
    fs_new = 500;
    
    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = 'selAreas_avgArea';
    savename_datestr_format = 'yyyymmdd';
    
    %% Code Start Here
    disp('Loading data_segments from Ying large file')
    load(fullfile(input_folder, input_file), 'data_segments');
    data_seg_origin = data_segments;
    clear data_segments


    %%% extract data_segments for each dateBK %%%
    disp('Extracing rest lfp for each DateBK')
    dateBKStr_prev = '';
    for di = 1 : length(data_seg_origin)
        data_1seg = data_seg_origin(di);

        if data_1seg.is_stn_dbs || data_1seg.is_gpi_dbs || ~strcmp(data_1seg.sleep_state, 'wake')
            clear data_1seg
            continue
        end

        dateBKStr = data_1seg.date;
        if ~strcmp(dateBKStr, dateBKStr_prev) % a new dateBK
            
            if ~strcmp(dateBKStr_prev, '') && ~isempty(data_segments)% save data using dateBKStr_prev
                tmps = split(dateBKStr_prev, '_');
                dateofexp = datenum(tmps{1}, 'yyyymmdd');
                bktdt = str2num(tmps{2});
                
                pdcond = parsePDCondition(dateofexp, animal);
                
                fs = fs_new;
                
                savefilename =  [animal '_' savefilename_addstr '_' pdcond '_' datestr(dateofexp, savename_datestr_format) '_tdt' num2str(bktdt) '.mat'];
                save(fullfile(savefolder, savefilename), 'data_segments', 'T_chnsarea', 'fs')
                
                clear tmps dateofexp bktdt pdcond savefilename
            end
            clear data_segments
            dateBKStr_prev = dateBKStr;

            % new data_segments and segi for the new dateBK
            data_segments = struct();
            segi = 1;
        end
        
        % extract lfp data: ntemp * nchns
        lfp_m1 = data_1seg.m1;
        lfp_stn = data_1seg.lfp_stn_diff;
        lfp_gp = data_1seg.lfp_gp_diff;
        
        % down sample 
        lfp_m1 = resample(lfp_m1, round(fs_new), round(fs_infile));
        lfp_stn = resample(lfp_stn, round(fs_new), round(fs_infile));
        lfp_gp = resample(lfp_gp, round(fs_new), round(fs_infile));
        
        % add chn-area information T_chnsarea
        if ~exist('T_chnsarea', 'var')
            nM1 = size(lfp_m1, 2);
            nSTN = size(lfp_stn, 2);
            nGP = size(lfp_gp, 2);
            T_chnsarea = chanInf_M1DBS(nM1, nSTN, nGP);

        else
            if nM1 ~= size(lfp_m1, 2) || nSTN ~= size(lfp_stn, 2) || nGP ~= size(lfp_gp, 2)
                disp(['nM1, nSTN and nGP not consistent for di = ' num2str(di)])
            end
            
        end

        % combine lfp_m1, lfp_stn and lfp_gp
        data_segments(segi).lfp = cat(1, lfp_m1', lfp_stn', lfp_gp');
        segi = segi + 1;

        clear data_1seg lfp_m1 lfp_stn lfp_gp 
    end
    
    % the last dateBKStr 
    if  ~isempty(data_segments)
        tmps = split(dateBKStr, '_');
        dateofexp = datenum(tmps{1}, 'yyyymmdd');
        bktdt = str2num(tmps{2});
        
        pdcond = parsePDCondition(dateofexp, animal);
        
        savefilename =  [animal '_' savefilename_addstr '_' pdcond '_' datestr(dateofexp, savename_datestr_format) '_tdt' num2str(bktdt) '.mat'];
        save(fullfile(savefolder, savefilename), 'data_segments', 'T_chnsarea', 'fs')
        
        clear tmps dateofexp bktdt pdcond savefilename
    end
    
    
    %%% remove unwanted chns %%% 
    disp('Removing unwanted chns')
    files = dir(fullfile(savefolder, '*.mat'));
    for filei = 1:length(files)
        load(fullfile(files(filei).folder, files(filei).name), 'data_segments', 'fs', 'T_chnsarea');
        
        
        % replace STN into stn0-1, stn1-2 ...
        row_STN = find(strcmp(T_chnsarea.brainarea, 'STN'));
        for i = 1: length(row_STN)
            recChn = T_chnsarea.recordingchn(row_STN(i));
            T_chnsarea.brainarea(row_STN(i)) = {['stn' num2str(recChn-1) '-' num2str(recChn)]};
            clear recChn
        end
        % replace GP into gp0-1, gp1-2 ...
        row_GP = find(strcmp(T_chnsarea.brainarea, 'GP'));
        for i = 1: length(row_GP)
            recChn = T_chnsarea.recordingchn(row_GP(i));
            T_chnsarea.brainarea(row_GP(i)) = {['gp' num2str(recChn-1) '-' num2str(recChn)]};
            clear recChn
        end
        
        % rm unwanted Chns in T_chnsarea and data_segments
        mask_remainChns = ~cellfun(@(x) any(strcmp(unwanted_DBS, x)), T_chnsarea.brainarea);
        T_chnsarea = T_chnsarea(mask_remainChns, :);
        T_chnsarea.chni = uint8([1:height(T_chnsarea)]');
        for segi = 1 : length(data_segments)
            lfp = data_segments(segi).lfp;
            data_segments(segi).lfp = lfp(mask_remainChns, :);
            clear lfp
        end
        
        
        % save to the same file
        save(fullfile(files(filei).folder, files(filei).name), 'data_segments', 'fs', 'T_chnsarea');

        clear data_segments fs T_chnsarea 
        clear row_STN row_GP mask_remainChns 
    end
    
 
end


function [T_chnsarea] = chanInf_M1DBS(nM1, nSTN, nGP)
    % add M1 and DBS channel together
    %
    %   Args:
    %       nM1: the number of M1 channels
    %       nstn, ngp: the number of stn and gp channels
    %
    %   Output:
    %       T_chnsarea: table of M1 + DBS channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    % initial the attitudes of chanarea table T_chnsarea:
    % chni_vec, electype, brainareas, notes, recordingchn
    chni_vec = uint8([1:nM1]');
    electype = cell(nM1, 1);
    brainareas = cell(nM1, 1);
    notes = cell(nM1, 1);
    recordingchn = uint8([1:nM1]');
    electype(1:nM1, 1) = {'Gray Matter'}; % electype
    brainareas(1:nM1, 1) = {'M1'}; % electype

    % channel information table of M1
    T_chnsarea = table;
    T_chnsarea.chni = chni_vec;
    T_chnsarea.brainarea = brainareas;
    T_chnsarea.recordingchn = recordingchn;
    T_chnsarea.electype = electype;
    T_chnsarea.notes = notes;

    % add DBS lead
    for chi = 1:nSTN
        T_chnsarea = [T_chnsarea; {nM1 + chi, 'STN', chi, 'DBS', ''}];
    end

    for chi = 1:nGP
        T_chnsarea = [T_chnsarea; {nM1 + nSTN + chi, 'GP', chi, 'DBS', ''}];
    end

end
    
    
    