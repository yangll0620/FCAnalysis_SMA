function m0_NormalCOTModerateSKTData_extract_Kitty()
    %% extract the SKT Data (normal COT and moderate STK data for Kitty)
    %
    %	Inputs:
    %
    %		xlsxfile_master
    %
    %   Steps:
    %       1. extract trials from processed folder in server (call function _extractlfptrial
    %           a. trial length = max(each trial length) + t_bef + t_aft 
    %                   t_bef: the time before target on (default: t_bef = 1)
    %                   t_aft: the time after mouth (default: t_aft = 0.5)
    %           b. only remain the trials markedin both good reach and good return
    %
    %
    %       2. bipolar for DBS channels
    %
    %       3. Down sample trials into fs_new = 500
    %
    %       4. add variable T_chnsarea


    %% extract the corresponding pipeline folder for this code
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');
    % code folder
    codefolder = codefilepath(1:strfind(codefilepath, 'code') + length('code') - 1);

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));
    addpath(genpath(fullfile(codefolder, 'NHPs')));

    % add NexMatablFiles path
    addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

    % datafolder, pipelinefolder
    [datafolder, ~, ~, ~] = exp_subfolders();
    correspipelinefolder = code_corresfolder(codefilepath, true, false);

    %% input setup
    animal = animal_extract(correspipelinefolder);

    % Input dir:  preprocessed folder in root2
    server_NHP = fullfile('Y:', 'root2', 'Animals');
    folder_datadatabase = fullfile(server_NHP, animal, 'Recording', 'Processed', 'DataDatabase');

    % master sheet
    xlsxfile_master = fullfile(datafolder, 'Animals', animal, 'KittyNormalCOTModerateSKT.csv');
    strformat_date_master = 'yyyymmdd'; % the format of date string in master sheet, e.g '20141013'


    % channel number of utah array, grayMatter, and dbs leads
    nM1 = 96; nSTN = 8; nGP = 8;


    % downsample
    fs_new = 500;
    
    % conds
    conds_cell = cond_cell_extract(animal);

    %% save setup
    savefolder = correspipelinefolder;
    savefilename_prefix = animal;

    strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'


    %% Starting Here

    % extract master table
    opts = detectImportOptions(xlsxfile_master);
    opts = setvartype(opts, {'Date'}, 'string');
    t_datebks = readtable(xlsxfile_master, opts);

    %%% add chn-area information T_chnsarea  %%%
    T_chnsarea = chanInf_M1DBS(nM1, nSTN -1, nGP -1); % bipolar DBS


    %%% extract all STK/COT trials %%%
    f = waitbar(0, 'Extracting all STK trials');
    nrecords = height(t_datebks);
    for i = 1:nrecords
        % waitbar
        waitbar(i / nrecords, f, ['Extracting trials in file ' num2str(i) '/' num2str(nrecords)]);

        % date of exp, bktdt
        dateofexp = datenum(t_datebks.Date(i), strformat_date_master);
        tdtbk = t_datebks.TDTBlock(i);
        

        % get the pd conditioon for the date of experiment
        pdcond = parsePDCondition(dateofexp, animal);
        

        if ~any(strcmp(pdcond, conds_cell)) % avoid tomild, tomoderate
            clear dateofexp tdtbk pdcond
            continue;
        end
        
        
        savefilename = [savefilename_prefix pdcond '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk) '.mat'];
        if exist(fullfile(savefolder, savefilename), 'file') == 2
            clear dateofexp tdtbk pdcond savefilename
            continue;
        end
         
        % one day path
        onedaypath = fullfile(folder_datadatabase, [animal '_' datestr(dateofexp, 'mmddyy')]);
        if ~exist(onedaypath, 'dir')
            clear dateofexp tdtbk pdcond savefilename onedaypath
            continue;
        end

        disp(['extracting ' onedaypath '-tdtbk' num2str(tdtbk)])
        
        % extract trials of lfp data for particular date and tdt block number
        if strcmp(pdcond, 'normal')
            [lfptrial_cortical, lfptrial_dbs, mask_goodreach, mask_goodreturn, fs_lfp, T_idxevent_lfp, fs_ma, T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial, T_dbsChn] = ...
                extractlfptrial_seg_normalCOT_Kitty(onedaypath, tdtbk);
        end
        if strcmp(pdcond, 'moderate')
            [lfptrial_cortical, lfptrial_dbs, mask_goodreach, mask_goodreturn, fs_lfp, T_idxevent_lfp, fs_ma, T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial, T_dbsChn] = ...
                extractlfptrial_seg(onedaypath, tdtbk);
        end
        
        
        % skip this day if any is empty
        if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs_lfp) || isempty(T_idxevent_lfp) || isempty(T_dbsChn)
            clear outputfoldername dateofexp tdtbk
            clear pdcond savefilename onedaypath
            clear lfptrial_cortical lfptrial_dbs fs_lfp T_idxevent_lfp mask_goodreach mask_goodreturn
            clear fs_ma T_idxevent_ma smoothWspeed_trial Wpos_smooth_trial Wrist_smooth_trial T_dbsChn
            
            continue;
        end
        
        
        %%% extract the lfptrial_m1 marked with good channels %%%
        lfptrial_cortical = lfptrial_cortical;
        
        %%% bipolar dbs %%%%
        for tri = 1 : length(lfptrial_dbs)
            lfp_1trial_dbs = lfptrial_dbs{tri};
            lfptrial_stn = diff(lfp_1trial_dbs(:, 1:nSTN), [], 2);
            lfptrial_gp = diff(lfp_1trial_dbs(:, nSTN + 1:end), [], 2);
            lfp_1trial_dbs = cat(2, lfptrial_stn, lfptrial_gp);
            lfptrial_dbs{tri} = lfp_1trial_dbs;
            clear lfp_1trial_dbs lfptrial_stn lfptrial_gp
        end
        clear tri
        
         % concatenate the utah chns, and the dbs chns into lfptrial
         lfpdata = cell(length(lfptrial_dbs), 1);
         for tri = 1 : length(lfptrial_dbs)
             lfpdata{tri} = cat(2, lfptrial_cortical{tri}, lfptrial_dbs{tri});
         end
        
        clear lfptrial_cortical lfptrial_dbs
        
        
        %%% down sample %%%
        for tri = 1 : length(lfpdata)
            lfpdata_1trial = lfpdata{tri};
            
            nchns = size(lfpdata_1trial, 2);
            for chi = 1 : nchns
                tmp = lfpdata_1trial(:, chi);
                tmp = resample(tmp, round(fs_new), round(fs_lfp));
                if chi == 1
                    lfpdown = zeros(nchns, length(tmp));
                end
                lfpdown(chi, :) = tmp;
            end
            lfpdata{tri} = lfpdown;
            clear lfpdata_1trial lfpdown nchns
        end
        T_idxevent_lfp{:, :} =  round(T_idxevent_lfp{:, :} * fs_new / fs_lfp);
        fs_lfp = fs_new;
        clear lfpdown nchns chi
        
        
         % save
        save(fullfile(savefolder, savefilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'T_idxevent_lfp', 'fs_ma', 'T_idxevent_ma', ...
            'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial', 'mask_goodreach', 'mask_goodreturn');
        
        
        clear outputfoldername dateofexp tdtbk
        clear pdcond savefilename onedaypath
        clear lfptrial_cortical lfptrial_dbs fs_lfp T_idxevent_lfp
        clear fs_ma T_idxevent_ma smoothWspeed_trial Wpos_smooth_trial Wrist_smooth_trial T_dbsChn
        clear lfpdata

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