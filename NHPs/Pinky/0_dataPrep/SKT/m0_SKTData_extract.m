function m0_SKTData_extract()
    %% extract the STK data  using spreetsheet
    %
    %	Inputs:
    %
    %		datafolder/Pinky/MasterDatabase.xlsx:  for stk file information
    %
    %   Steps:
    %       1. extract trials from processed folder in server (call function _extractlfptrial
    %           a. trial length = max(each trial length) + t_bef + t_aft 
    %                   t_bef: the time before target on (default: t_bef = 1)
    %                   t_aft: the time after mouth (default: t_aft = 0.5)
    %           b. only remain the trials markedin both good reach and good return
    %
    %       2. only keep the good M1 channels (chans_m1_remain.mat) and GM channels marked with brain area
    %
    %       3. bipolar for DBS channels
    %
    %       4. Down sample trials into fs_new = 500
    %
    %       5. add variable T_chnsarea


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
    animal = 'Pinky';

    % Input dir:  preprocessed folder in root2
    processedfolder_inroot2 = fullfile('H:', 'My Drive', 'NMRC_umn', 'Projects', 'FCAnalysis', 'exp', 'data',animal, 'Recording', 'Processed', 'DataDatabase');

    % master sheet
    xlsxfile_master = fullfile(datafolder, animal, [animal 'MasterDatabase.xlsx']);
    strformat_date_master = 'mmddyy'; % the format of date string in master sheet, e.g '012317'

    % different SKB labels in the master sheet
    tasks_SKB = {'SKB', 'Single', 'Single KB', 'Single Kluver', 'Single Target', 'Single Target Kluver', 'Single-Target Kluver', ...
                'Single-target Kluver', 'single', 'single Kluver', 'single-target Kluver'};

    % good M1 channels stored file
    matfile_m1chns = fullfile(datafolder, animal, 'chans_m1_remain.mat');

    % gray matter chn-area information file
    filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
    file_GMChnsarea = fullfile(datafolder, animal, filename_GMChnsarea);

    % channel number of utah array, grayMatter, and dbs leads
    nM1 = 96; nGM = 32; nSTN = 8; nGP = 8;


    % downsample
    fs_new = 500;
    
    % conds
    conds_cell = {'normal', 'mild', 'moderate'};

    %% save setup
    savefolder = correspipelinefolder;
    savefilename_prefix = [animal '_STKData_'];

    strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'


    %% Starting Here
    % load good channels of M1 (i.e. Utah array)
    load(matfile_m1chns, 'chans_m1');

    % extract master table
    t_master = readtable(xlsxfile_master);

    % find the row indices for task_SKB, idx_rows: 0 for not task_SKB, 1 for task_SKB
    idx_rows = cellfun(@(x) ~isempty(find(strcmp(x, tasks_SKB))), t_master.BriefDescription);

    % table for rows marked SKB labels with only OutputFolderName and TDTBlock columns
    t_SKT = t_master(idx_rows, {'OutputFolderName', 'TDTBlock'});

    % convert Table Variables from Cell Arrays of Character Vectors to string/double Arrays
    t_SKT.OutputFolderName = string(t_SKT.OutputFolderName);
    t_SKT.TDTBlock = double(t_SKT.TDTBlock);

    %%% add chn-area information T_chnsarea  %%%
    T_chnsarea_M1 = chanInf_M1(chans_m1); % T_chnsarea_M1

    T_chnsarea_GM = chanInf_GM(file_GMChnsarea, nGM); % GM chn Areas table
    % only keep the GM channels with brainarea
    chns_GM = cellfun(@(x)~isempty(x), T_chnsarea_GM.brainarea);
    T_chnsarea_GM = T_chnsarea_GM(chns_GM, :);

    T_chnsarea_DBS = chanInf_DBS(nSTN - 1, nGP -1); % bipolar DBS

    % adjust chni column
    T_chnsarea_GM.chni = T_chnsarea_GM.chni + nM1;
    T_chnsarea_DBS.chni = T_chnsarea_DBS.chni + nM1 + nGM;
        
    T_chnsarea = [T_chnsarea_M1; T_chnsarea_GM; T_chnsarea_DBS];
    clear T_chnsarea_M1 T_chnsarea_DBS T_chnsarea_GM


    %%% extract all STK trials %%%
    f = waitbar(0, 'Extracting all STK trials');
    n = height(t_SKT);

    for i = 1:n
        % waitbar
        waitbar(i / n, f, ['Extracting trials in file ' num2str(i) '/' num2str(n)]);

        % date of exp, bktdt
        outputfoldername = split(t_SKT.OutputFolderName(i), '_');
        dateofexp = datenum(outputfoldername{end}, strformat_date_master);
        tdtbk = t_SKT.TDTBlock(i);


        % get the pd conditioon for the date of experiment
        pdcondition = parsePDCondition(dateofexp, animal);
        
        if ~strcmp(pdcondition, 'moderate')
            continue;
        end

        if ~any(strcmp(pdcondition, conds_cell)) % avoid tomild, tomoderate
            continue;
        end
         
        % one day path
        onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);

        disp(['extracting ' onedaypath '-tdtbk' num2str(tdtbk)])

        % extract trials of lfp data for particular date and tdt block number
        [lfptrial_cortical, lfptrial_dbs, fs, T_idxevent, T_dbsChn] = extractlfptrial_(onedaypath, tdtbk);

        % skip this day if any is empty
        if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) || isempty(T_dbsChn)
            continue;
        end
        

        %%% extract the lfptrial_m1 marked with good channels %%%
        lfptrial_m1 = lfptrial_cortical(1:nM1, :, :);
        lfptrial_m1 = lfptrial_m1(chans_m1, :, :);

        %%% extract the lfptrial_GM used in brain area %%%
        lfptrial_GM = lfptrial_cortical(end - nGM + 1:end, :, :);
        lfptrial_GM = lfptrial_GM(chns_GM, :, :);

        %%% bipolar dbs %%%%
        lfptrial_stn = diff(lfptrial_dbs(1:nSTN, :, :), [], 1);
        lfptrial_gp = diff(lfptrial_dbs(nSTN+1:nSTN+nGP, :, :), [], 1);
        lfptrial_dbs = cat(1, lfptrial_stn, lfptrial_gp);


        %%% down sample %%%
        nchns = size(lfptrial_m1, 1);
        for chi = 1: nchns
            tmp = squeeze(lfptrial_m1(chi, :, :));
            tmp = resample(tmp, round(fs_new), round(fs));
            
            if chi == 1
                [s1, s2] = size(tmp);
                lfpdown_m1 = zeros(nchns, s1, s2);
                clear s1 s2
            end
            lfpdown_m1(chi, :, :) = tmp;
             
            clear tmp
        end
        nchns = size(lfptrial_GM, 1);
        for chi = 1: nchns
            tmp = squeeze(lfptrial_GM(chi, :, :));
            tmp = resample(tmp, round(fs_new), round(fs));
            
            if chi == 1
                [s1, s2] = size(tmp);
                lfpdown_GM = zeros(nchns, s1, s2);
                clear s1 s2
            end
            lfpdown_GM(chi, :, :) = tmp;
             
            clear tmp
        end
        nchns = size(lfptrial_dbs, 1);
        for chi = 1: nchns
            tmp = squeeze(lfptrial_dbs(chi, :, :));
            tmp = resample(tmp, round(fs_new), round(fs));
            
            if chi == 1
                [s1, s2] = size(tmp);
                lfpdown_dbs = zeros(nchns, s1, s2);
                
                clear s1 s2
            end
            lfpdown_dbs(chi, :, :) = tmp;
            
            clear tmp
        end
        T_idxevent{:, :} =  round(T_idxevent{:, :} * fs_new / fs);
        fs = fs_new;


        % concatenate the utah chns, GM chns and the dbs chns into lfptrial
        lfpdata = cat(1, lfpdown_m1, lfpdown_GM, lfpdown_dbs);


        % save
        savefilename = [savefilename_prefix pdcondition '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk)];
        save(fullfile(savefolder, savefilename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');

        clear outputfoldername dateofexp tdtbk onedaypath
        clear lfptrial_cortical lfptrial_dbs fs T_idxevent T_dbsChn lfptrial_m1 lfptrial_GM
        clear lfpdata  pdcondition savefilename
    end
    
    close(f)
end

function [lfptrial_cortical, lfptrial_dbs, fs, idxeventtbl, chantbl_dbs] = extractlfptrial_(onedaypath, tdtblock)
    % extractlfptrial extract trials for LFP data
    %
    %  [lfptrial_cortical, lfptrial_dbs, chantbl_cortical, chantbl_dbs] =
    %  extractlfptrial(onedaypath, block) return extracted LFP trials of
    %  cortical/subcortical data, dbs data, cortical/subcortical channel
    %  information and dbs channel information tables, only trials with both
    %  good reach and return are returned
    %
    %  Example usage:
    %   onedaypath = '/home/lingling/root2/Animals2/Pinky/Recording/Processed/DataDatabase/Pinky_091517'
    %   tdtblock = 3
    %   [lfptrial_cortical, lfptrial_dbs, idxtbl_event,chantbl_cortical, chantbl_dbs] = extractlfptrial(onedaypath, tdtblock)
    %
    %  Inputs:
    %   onedaypath : one date folder
    %   tdtblock: tdt block number
    %
    %  Used files:
    %       lfpfile_utah  -   .\LFP\Block-3\Pinky_GrayMatter_eyetracking_DT1_091517_Block-3_LFPch*.nex
    %       lfpfile_dbs   -   .\DBSLFP\Block-3\Pinky_GrayMatter_eyetracking_DT1_091517_Block-3_DBSLFP.nex
    %       mafile        -   .\Block-3\pinky_20170915_4_cleaned_MA_SingleTargetKluver_Analyze2.mat
    %
    %  Outputs:
    %   lfptrial_cortical: lfp trials of cortical/subcortical channels
    %                      [chn_cortical * n_temporal * n_trial]
    %
    %        lfptrial_dbs: lfp trials of dbs channels
    %                      [chn_dbs * n_temporal * n_trial], 1-8: STN, 9-16:GP
    %        fs: sample rate
    %
    %
    %        idxeventtbl: a table describes the index for events of target onset,
    %                      reach onset, touch screen, return and mouth in the trial
    %
    %         chantbl_dbs:  a table describes each dbs channel information
    %
    %  More Description:
    %       trial length = max(each trial length) + t_bef + t_aft
    %       t_bef: the time before target on (default: t_bef = 1)
    %       t_aft: the time after mouth (default: t_aft = 0.5)
    %       one trailis from 't_target - t_bef'  to 't_mouth + t_after'

    %% add NexMatablFiles path

    [~, codefolder, ~, ~] = exp_subfolders();
    addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

    %% MA data
    % read the MA data
    mafolder = fullfile(onedaypath, ['Block-' num2str(tdtblock)]); 
    mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));

    % ma file does not exist
    if length(mafilestruct) ~= 1
        disp([mafolder ' has ' num2str(length(mafilestruct)) ' files, skip!'])

        lfptrial_cortical = []; lfptrial_dbs = [];
        fs = [];
        idxeventtbl = []; chantbl_dbs = [];

        return;
    end

    % load SingleTargetKluverMAData
    load(fullfile(mafolder, mafilestruct.name), 'SingleTargetKluverMAData');

    % ma sample rate
    fs_ma = SingleTargetKluverMAData.SR;

    % time indices for target onset, reach onset, touch screen, return and mouth
    TargetTime = SingleTargetKluverMAData.TargetTime;
    ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
    TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
    ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
    MouthTimeix = SingleTargetKluverMAData.MouthTimeix;

    timeixtbl_ma = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];

    % the tag of good reach trials
    tag_goodreach = SingleTargetKluverMAData.goodix_reach;
    % the tag of good return trials
    tag_goodreturn = SingleTargetKluverMAData.goodix_return;

    % extract indices of good trials (both have good reach and return)
    idx_goodtrials = (tag_goodreach == 1) & (tag_goodreturn == 1);

    % if no good trials can be found
    if isempty(idx_goodtrials)
        disp(['no good trials are found in ' mafilestruct.name])

        lfptrial_cortical = []; lfptrial_dbs = [];
        fs = [];
        idxeventtbl = []; chantbl_dbs = [];

        return;
    end

    %
    timeixtbl_ma = timeixtbl_ma(idx_goodtrials, :);

    clear TargetTime ReachTimeix TouchTimeix ReturnTimeix MouthTimeix

    %% LFP data
    t_bef = 1; t_aft = 0.5; % t_bef: time before target on, t_aft: time after mouth
    [n_trial, n_timevars] = size(timeixtbl_ma);

    % read each channel data in  LFP data to lfpdata (initialized when is the 1st channel)
    folder_cortical = fullfile(onedaypath, 'LFP', ['Block-' num2str(tdtblock)]);

    % files in folder_cortical are empty
    if isempty(dir(folder_cortical))
        disp([folder_cortical ' has no files!'])

        lfptrial_cortical = []; lfptrial_dbs = [];
        fs = [];
        idxeventtbl = []; chantbl_dbs = [];

        return;
    end

    filenames = extractfield(dir(folder_cortical), 'name');
    match = cellfun(@(x)~isempty(regexp(x, ['LFPch[0-9]*.nex'], 'match')), filenames, 'UniformOutput', false); % match channel file
    match = cell2mat(match);
    nexnames = filenames(match);

    % parse and sort the number of channel
    for filei = 1:length(nexnames)
        filename = nexnames{filei};
        tmp = char(regexp(filename, 'ch[0-9]*+.nex', 'match'));
        chns(filei) = str2num(tmp(length('ch') + 1:end - length('.nex')));

        if filei == 1
            % get the prefix of the each nex file name
            tmp = split(filename, [num2str(chns(filei)) '.nex']);
            file_prefix = tmp{1};
        end

    end

    chns = sort(chns);
    chn_lfp = length(chns);

    for i = 1:length(chns)
        filename = [file_prefix num2str(chns(i)) '.nex'];
        [nexlfp_cortical] = readNexFile(fullfile(folder_cortical, filename));

        % extract the number of the structure containing LFPchn* data
        name_list = extractfield(cell2mat(nexlfp_cortical.contvars), 'name');
        i_lfp = find(contains(name_list, 'LFP')); %% i.e nexlfp_utah.contvars name == 'LFPch1', or 'MUAch1'

        if i == 1% first channel
            fs_lfpcortical = nexlfp_cortical.contvars{i_lfp}.ADFrequency;

            % the time index in the LFP neural data based on MA data
            timeixtbl_lfpcortical = timeixtbl_ma;
            timeixtbl_lfpcortical{:, :} = round(timeixtbl_ma{:, :} / fs_ma * fs_lfpcortical);

            % initialize lfp_utah : chn_lfp * n_temporal * n_trial
            maxlen = max(timeixtbl_lfpcortical.MouthTimeix - timeixtbl_lfpcortical.TargetTime) + 1; % maximum length across all trials (unit: ind)
            n_bef = round(t_bef * fs_lfpcortical); % n_bef: index number before target on
            n_aft = round(t_aft * fs_lfpcortical); % n_aft: index number after mouth
            n_temporal = maxlen + n_bef + n_aft;
            idx_str = timeixtbl_lfpcortical.TargetTime - n_bef;
            idx_end = timeixtbl_lfpcortical.TargetTime + (maxlen -1) + n_aft;
            lfptrial_cortical = zeros(chn_lfp, n_temporal, n_trial);

            % idxtbl_lfptrialutah: the idx for events of target onset, reach onset, touch screen,
            % return and mouth in the trial matrix lfpdata_utah (chn_lfputah * n_temporal * n_trial)
            % the first sample corresponds to target onset - t_bef
            idxtbl_lfptrial_cortical = timeixtbl_lfpcortical;
            idxtbl_lfptrial_cortical{:, :} = idxtbl_lfptrial_cortical{:, :} - repmat(idx_str, [1, n_timevars]);

            clear timeixtbl_lfpsepchn n_bef n_aft
        else

            if fs_lfpcortical ~= nexlfp_cortical.contvars{i_lfp}.ADFrequency% samping frequency is different
                chni = chns(i);
                disp(['sampling frequency is different for chni = ' num2str(chni)]);
                break;
            end

        end

        % extract each trial for lfp data stored in separate channel
        for triali = 1:n_trial

            if size(nexlfp_cortical.contvars{i_lfp}.data, 1) < idx_end(triali)
                disp(mafilestruct)
                disp(filename)
                disp(fs_ma)
                disp(fs_lfpcortical)
                disp(['idx in ma' num2str(timeixtbl_ma{1, end})])
                disp(size(nexlfp_cortical.contvars{i_lfp}.data, 1))
                disp(idx_end(triali))
                disp(['maxlen  = ' num2str(maxlen)])
            end

            lfptrial_cortical(i, :, triali) = nexlfp_cortical.contvars{i_lfp}.data(idx_str(triali):idx_end(triali));

        end

        clear filename i_lfp
    end

    % disp play the max trial time
    disp(['max trial time is ' num2str(maxlen / fs_lfpcortical)]);


    %% DBSLFP data
    % read DBSLFP in nex file
    dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtblock)]); % 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\DBSLFP\Block-1'

    % dbs file does not exist
    if isempty(dir(fullfile(dbslfpfolder, '*_GrayMatter*DBSLFP.nex')))
        disp([dbslfpfolder, '*_GrayMatter*DBSLFP.nex files do not exit!'])

        lfptrial_cortical = []; lfptrial_dbs = [];
        fs = [];
        idxeventtbl = []; chantbl_dbs = [];

        return;
    end

    filenames = extractfield(dir(fullfile(dbslfpfolder, '*_GrayMatter*DBSLFP.nex')), 'name'); % Pinky_GrayMatter_eyetracking_DT1_071417_Block-1_DBSLFP.nex
    match = cellfun(@(x)~isempty(regexp(x, ['[A-Z, a-z, 0-9]*.nex'], 'match')), filenames, 'UniformOutput', false); % match channel file
    match = cell2mat(match);
    nexnames = filenames(match);

    if length(nexnames) ~= 1
        % only has one file inside folder
        error(['More than one _DBSLFP.nex  in file' dbslfpfolder])
    end

    [nexlfp_dbs] = readNexFile(fullfile(dbslfpfolder, nexnames{1}));

    % find the STN and GP data in the struct nexlfp_dbs
    % extract all the values of field 'name' in nexlfp_dbs.convars
    convars = cell2mat(nexlfp_dbs.contvars);
    name_list = extractfield(cell2mat(nexlfp_dbs.contvars), 'name');
    % extract the indices representing STN and GP data in nexlfp_dbs.contvars (33*1 struct array)
    idx_dbs = find(contains(name_list, 'RAW_DBSch'));

    % check the dbs channel frequencies
    if range(cell2mat({convars(idx_dbs).ADFrequency})) ~= 0
        % check the sampling frequencies in channels are consistent
        disp(['ADFrequency in the dbs channels (STN & GP) is not consistent']);
        return
    end

    if (convars(idx_dbs(1)).ADFrequency ~= fs_lfpcortical)
        % check the sampling frequency of dbs is the same as the lfp stored in separate channels or not
        disp(['the sampling frequency of dbs is not the same as the lfp stored in separate channels'])
        return
    end

    fs = fs_lfpcortical; % sampling frequency for dbs channels
    idxeventtbl = idxtbl_lfptrial_cortical;
    clear idxtbl_lfptrial_sepchn

    % extract each trial for lfp dbs data
    nchn_dbs = length(idx_dbs);
    lfptrial_dbs = zeros(nchn_dbs, n_temporal, n_trial);

    for i = 1:nchn_dbs
        chni = idx_dbs(i);
        lfp_1chn = convars(chni).data;

        for triali = 1:n_trial
            lfptrial_dbs(i, :, triali) = lfp_1chn(idx_str(triali):idx_end(triali));
        end

        clear chni lfp_1chn triali
    end

    % chantbl_dbs: a table describes each dbs channel information
    elecchn = extractfield(convars(idx_dbs), 'name');
    elecchn = elecchn';
    area = cell(nchn_dbs, 1);
    area(1:8) = {'STN'};
    area(9:16) = {'GP'};
    chantbl_dbs = [table(area) table(elecchn)];
    clear varName

    clear convars filename idx_stn idx_gp
end

function T_chnsarea = chanInf_GM(file_GMChnsarea, nGM)
    % extract gray matter channel inf table (recording thalamus, SMA et al. areas)
    %
    %   Args:
    %       file_GMCchnsarea: the file storing the gray matter chn-area inf
    %       nGM: the total number of gray matter channels
    %
    %   Output:
    %       T_chnsarea: table of Gray matter channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    T = readtable(file_GMChnsarea);
    chi_firstGM = 101;
    nareas = height(T);
    GMChnAreas = cell(nareas, 1);

    for areai = 1:nareas
        area = T.brainarea{areai};
        tmpcell = split(T.channels{areai}, ',');

        for j = 1:length(tmpcell)
            chn = str2num(char(tmpcell{j}));
            GMChnAreas(chn - chi_firstGM + 1, 1) = {area};
        end

    end

    % channel information table
    T_chnsarea = table;
    T_chnsarea.chni = uint8([1:nGM]');
    T_chnsarea.brainarea = GMChnAreas;
    T_chnsarea.recordingchn = uint8([1:nGM]') + chi_firstGM -1;
    T_chnsarea.electype = cell(nGM, 1);
    T_chnsarea.electype(:) = {'Gray Matter'};
    T_chnsarea.notes = cell(nGM, 1);

end

function T_chnsarea = chanInf_M1(chans_m1)
    % extract M1 channel inf table
    %   Args:
    %       chans_m1: a vector containing m1 good channel numbers
    %
    %   Output:
    %       T_chnsarea: table of M1 channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    nM1 = length(chans_m1);
    chans_m1 = reshape(chans_m1, nM1, 1);

    % channel information table of M1
    T_chnsarea = table;
    T_chnsarea.chni = uint8([1:nM1]');
    T_chnsarea.brainarea = cell(nM1, 1);
    T_chnsarea.brainarea(:) = {'M1'};
    T_chnsarea.recordingchn = chans_m1;
    T_chnsarea.electype = cell(nM1, 1);
    T_chnsarea.electype(:) = {'Utah Array'};
    T_chnsarea.notes = cell(nM1, 1);

end

function T_chnsarea = chanInf_DBS(nSTN, nGP)
    % extract M1 channel inf table
    %   Args:
    %       nSTN, nGP: the number of stn and gp channels
    %
    %   Output:
    %       T_chnsarea: table of DBS channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    % channel information table of M1
    T_chnsarea = table;

    T_chnsarea.chni = uint8([1:nSTN + nGP]');

    T_chnsarea.brainarea = cell(nSTN + nGP, 1);
    T_chnsarea.brainarea(1:nSTN) = {'STN'}; T_chnsarea.brainarea(nSTN + 1:nSTN + nGP) = {'GP'};

    T_chnsarea.recordingchn = [uint8([1:nSTN]'); uint8([1:nGP]')];

    T_chnsarea.electype = cell(nSTN + nGP, 1);
    T_chnsarea.electype(:) = {'DBS'};

    T_chnsarea.notes = cell(nSTN + nGP, 1);

end
