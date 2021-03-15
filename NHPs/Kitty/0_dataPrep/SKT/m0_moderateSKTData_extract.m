function m0_moderateSKTData_extract()
    %% extract the SKT moderate data for Kitty 
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

    % master sheet: Kitty doesn't have SKT in normal
    SKTfile_master = fullfile(datafolder, animal, 'KittyModerateSKT.csv');
    strformat_date_master = 'yyyymmdd'; % the format of date string in  SKTfile_master, e.g. '20170123'


    % channel number of utah array, grayMatter, and dbs leads
    nM1 = 96; nSTN = 8; nGP = 8;


    % downsample
    fs_new = 500;
    
    conds_cell = {'moderate'};
    
    strformat_date_onedaypath = 'mmddyy'; % the format of date string in  savedaypath, e.g. '012317'

    %% save setup
    savefolder = correspipelinefolder;
    savefilename_prefix = [animal '_STKData_'];

    strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'


    %% Starting Here

    % extract master table
    t_SKT = readtable(SKTfile_master, 'Format','%s%u%u%s');


    %%% add chn-area information T_chnsarea  %%%
    T_chnsarea = chanInf_M1DBS(nM1, nSTN -1, nGP -1); % bipolar DBS


    %%% extract all STK trials %%%
    f = waitbar(0, 'Extracting all STK trials');
    nrecords = height(t_SKT);

    for i = 1:nrecords
        % waitbar
        waitbar(i / nrecords, f, ['Extracting trials in file ' num2str(i) '/' num2str(nrecords)]);

        % date of exp, bktdt, condition
        dateofexp = datenum(t_SKT.Date(i), strformat_date_master);
        tdtbk = t_SKT.TDTBlock(i);
        pdcondition = t_SKT.Condition{i};
        

        if ~any(strcmp(pdcondition, conds_cell)) % avoid tomild, tomoderate
            continue;
        end
        
        savefilename = [savefilename_prefix pdcondition '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk) '.mat'];
        if exist(fullfile(savefolder, savefilename), 'file') == 2
            continue;
        end
         
        % one day path
        onedaypath = fullfile(processedfolder_inroot, [animal '_' datestr(dateofexp, strformat_date_onedaypath)]);

        disp(['extracting ' onedaypath '-tdtbk' num2str(tdtbk)])

        % extract trials of lfp data for particular date and tdt block number
        [lfptrial_cortical, lfptrial_dbs, fs, T_idxevent, T_dbsChn] = extractlfptrial_(onedaypath, tdtbk);

        % skip this day if any is empty
        if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(T_idxevent) || isempty(T_dbsChn)
            continue;
        end
        

        %%% extract the lfptrial_m1 marked with good channels %%%
        lfptrial_m1 = lfptrial_cortical;

        
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
            end
            lfpdown_m1(chi, :, :) = tmp;
             
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


        % concatenate the utah chns, and the dbs chns into lfptrial
        lfpdata = cat(1, lfpdown_m1, lfpdown_dbs);


        % save
        savefilename = [savefilename_prefix pdcondition '_' datestr(dateofexp, strformat_date_save) '_bktdt' num2str(tdtbk) '.mat'];
        save(fullfile(savefolder, savefilename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');

        clear outputfoldername dateofexp tdtbk onedaypath
        clear lfptrial_cortical lfptrial_dbs fs T_idxevent T_dbsChn lfptrial_m1 
        clear lfpdata  pdcondition savefilename
    end

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
    %   onedaypath = '/home/lingling/root2/Animals/Jo/Recording/Processed/DataDatabase/Jo_020416'
    %   tdtblock = 3
    %   [lfptrial_cortical, lfptrial_dbs, fs, idxeventtbl, chantbl_dbs] = extractlfptrial_(onedaypath, tdtblock);
    %
    %  Inputs:
    %   onedaypath : one date folder
    %   tdtblock: tdt block number
    %
    %  Used files:
    %       lfpfile_utah  -   .\LFP\Block-3\Jo_CR1_DT1_020216_Block-3_LFPch*.nex
    %       lfpfile_dbs   -   .\DBSLFP\Block-3\Jo_CR1_DT1_020216_Block-3_DBSLFP.nex.nex
    %       mafile        -   .\Block-3\Jo_20160202_3_cleaned_MA_SingleTargetKluver_Analyze2.mat
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
    %       one trailis from 't_target - t_bef'  to 't_mouth - t_after'

    %% add NexMatablFiles path

    [~, codefolder, ~, ~] = exp_subfolders();
    addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

    %% MA data
    % read the MA data
    mafolder = fullfile(onedaypath, ['Block-' num2str(tdtblock)]); 
    mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));

    % ma file does not exist
    if length(mafilestruct) ~= 1
        disp([mafolder 'has ' num2str(length(mafilestruct)) ' files, skip!'])

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
    TargetTime = SingleTargetKluverMAData.TargetTimeix;
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

%     if maxlen / fs_lfpcortical > ts_maxtrial
%         disp(['Abandon the file with max trial time > ' num2str(ts_maxtrial) 's'])
%         lfptrial_cortical = []; lfptrial_dbs = []; fs = []; idxeventtbl = []; chantbl_dbs = [];
%         return
%     end

    %% DBSLFP data
    % read DBSLFP in nex file
    dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtblock)]); % 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\DBSLFP\Block-1'

    dbsfilepattern = fullfile(dbslfpfolder, '*DBSLFP.nex');
    dbsfiles = dir(dbsfilepattern);
    
    % dbs file does not exist
    if length(dbsfiles) ~= 1
        disp([dbsfilepattern ' has ' num2str(length(dbsfiles)) ' file, skip!'])

        lfptrial_cortical = []; lfptrial_dbs = [];
        fs = [];
        idxeventtbl = []; chantbl_dbs = [];

        return;
    end

    [nexlfp_dbs] = readNexFile(fullfile(dbslfpfolder, dbsfiles(1).name));

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
