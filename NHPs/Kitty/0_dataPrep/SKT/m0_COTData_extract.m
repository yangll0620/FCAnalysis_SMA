function m0_COTData_extract()
    %% extract the COT normal data for Kitty 
    %
    %	Inputs:
    %
    %		root/Animals/Kitty/Recording/Processed/DataDatabase/
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

    % add NexMatablFiles path
    addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

    % datafolder, pipelinefolder
    [datafolder, ~, ~, ~] = exp_subfolders();
    correspipelinefolder = code_corresfolder(codefilepath, true, false);

    %% input setup
    tmp = regexp(codefilepath, 'NHPs/[a-zA-Z]*', 'match');
    animal = tmp{1}(length('NHPs/')+1:end);
    clear tmp

    % Input dir:  preprocessed folder in root2
    processedfolder_inroot = fullfile('/home', 'lingling', 'root', 'Animals', animal, 'Recording', 'Processed', 'DataDatabase');
    
    % SPK folder
    folder_SPK = fullfile(datafolder, animal, 'SPK');

    % master sheet: Kitty only has COT in normal 
    COTfile_master = fullfile(datafolder, animal, 'KittyNormalCOT.csv');
   
    % channel number of utah array, grayMatter, and dbs leads
    nM1 = 96; nSTN = 8; nGP = 8;


    % downsample
    fs_new = 500;
    
    conds_cell = {'normal'};
    
    strformat_date_master = 'yyyymmdd'; % the format of date string in  COTfile_master, e.g. '20170123'
    strformat_date_onedaypath = 'mmddyy'; % the format of date string in  savedaypath, e.g. '012317'
    strformat_date_SPKfilename = 'mmddyy'; % the format of date string in SPKfilename, e.g. '101614'

    %% save setup
    savefolder = correspipelinefolder;
    savefilename_prefix = [animal '_COTData_'];

    strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'






    %% Starting Here
    
    % extract master table
    t_COT = readtable(COTfile_master, 'Format','%s%u%u%s');
    
    
    
    %%% add chn-area information T_chnsarea  %%%
    T_chnsarea = chanInf_M1DBS(nM1, nSTN -1, nGP -1); % bipolar DBS
    
    
    %%% extract all COT trials %%%
    f = waitbar(0, 'Extracting all COT trials');
    nrecords = height(t_COT);
    for ni = 1:nrecords
        
        % waitbar
        waitbar(ni / nrecords, f, ['Extracting trials in file ' num2str(ni) '/' num2str(nrecords)]);

        
        % date of exp, bktdt, condition
        dateofexp = datenum(t_COT.Date(ni), strformat_date_master);
        tdtbk = t_COT.TDTBlock(ni);
        pdcondition = t_COT.Condition{ni};


        onedaypath = fullfile(processedfolder_inroot, [animal '_' datestr(dateofexp, strformat_date_onedaypath)]);

        
        % extract file_SPK containing MA information
        filenamepatt_SPK = ['*_' datestr(dateofexp, strformat_date_SPKfilename) '_Block-' num2str(tdtbk) '*.nex'];
        files_SPK = dir(fullfile(folder_SPK, filenamepatt_SPK));

        % SPK nex file containing MA information does not exist
        if length(files_SPK) ~= 1
            disp([folder_SPK 'has ' num2str(length(files_SPK)) ' files, skip!'])

            continue;
        end
        file_SPK = fullfile(folder_SPK, files_SPK.name);
        
        
        %extract trials of lfp data for particular date and tdt block number
        [lfptrial_notDBS, lfptrial_dbs, fs, idxeventtbl, chantbl_dbs] = extractlfptrial_COT(onedaypath, tdtblock, file_SPK);
        
        
        
        % skip this day if any is empty
        if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) || isempty(T_dbsChn)
            continue;
        end
        
        
        
        %%% bipolar dbs %%%%
        lfptrial_stn = diff(lfptrial_dbs(1:nSTN, :, :), [], 1);
        lfptrial_gp = diff(lfptrial_dbs(nSTN+1:nSTN+nGP, :, :), [], 1);
        lfptrial_dbs = cat(1, lfptrial_stn, lfptrial_gp);
        
        
        %%% cat notDBS and DBS LFP %%%
        lfptrial = cat(1, lfptrial_notDBS, lfptrial_dbs);
        
        %%% down sample %%%
        nchns = size(lfptrial, 1);
        for chi = 1: nchns
            tmp = squeeze(lfptrial(chi, :, :));
            tmp = resample(tmp, round(fs_new), round(fs));
            
            if chi == 1
                [s1, s2] = size(tmp);
                lfpdown = zeros(nchns, s1, s2);
            end
            lfpdown(chi, :, :) = tmp;
             
            clear tmp
        end
        T_idxevent{:, :} =  round(T_idxevent{:, :} * fs_new / fs);
        fs = fs_new;
        
        
        %%% final lfpdata %%%
        lfpdata = lfpdown;
        
        %%% save %%%
        savefilename = [savefilename_prefix pdcondition '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk) '.mat'];
        save(fullfile(savefolder, savefilename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');

        clear outputfoldername dateofexp tdtbk onedaypath
        clear lfptrial_cortical lfptrial_dbs fs T_idxevent T_dbsChn lfptrial_m1 
        clear lfpdata  pdcondition savefilename

    end

