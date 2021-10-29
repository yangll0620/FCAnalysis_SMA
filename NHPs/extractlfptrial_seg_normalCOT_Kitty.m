function [lfptrial_cortical, lfptrial_dbs, mask_goodreach, mask_goodreturn, fs_lfp, T_idxevent_lfp, fs_ma, T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial, T_dbsChn] = ...
    extractlfptrial_seg_normalCOT_Kitty(onedaypath, tdtbk)
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
    %       lfpfile_utah  -   .\LFP\Block-3\KittyArrayDBSv2_DT1_101614_Block-3_LFPch*.nex
    %       lfpfile_dbs   -   .\DBSLFP\Block-3\KittyArrayDBSv2_DT1_101614_Block-3_DBSLFP.nex
    %       mafile        -   .\Block-3\Kitty20141016_1_clean_MA_COT_Analyze2.mat
    %
    %  Outputs:
    %   lfptrial_cortical: lfp trials of cortical/subcortical channels
    %                      [chn_cortical * n_temporal * n_trial]
    %
    %        lfptrial_dbs: lfp trials of dbs channels
    %                      [chn_dbs * n_temporal * n_trial], 1-8: STN, 9-16:GP
    %        fs_lfp: sample rate for lfp data
    %
    %
    %        T_idxevent_lfp: a table describes the index for events of target onset,
    %                      reach onset, touch screen, return and mouth in the lfp trial
    %
    %        fs_ma: sample rate for ma data
    %
    %
    %        T_idxevent_ma: a table describes the index for events of target onset,
    %                      reach onset, touch screen, return and mouth in the ma data trial
    %
    %
    %        smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial: MA trial data
    %
    %        T_dbsChn:  a table describes each dbs channel information
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
    mafolder = fullfile(onedaypath, ['Block-' num2str(tdtbk)]);
    mafilestruct = dir(fullfile(mafolder, '*COT_Analyze2.mat'));

    if length(mafilestruct) ~= 1
        
        disp([mafolder 'has ' num2str(length(mafilestruct)) ' files, skip!'])
        
        lfptrial_cortical = []; lfptrial_dbs = [];
        mask_goodreach = []; mask_goodreturn = [];
        fs_lfp = []; T_idxevent_lfp = [];
        fs_ma = []; T_idxevent_ma = [];
        smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
        T_dbsChn = [];
        
        return
    end
    
    % load COTData from  *SingleTargetKluver_Analyze2.mat
    load(fullfile(mafolder, mafilestruct.name), 'COTData');

    % ma sample rate
    fs_ma = COTData.SR;

    % time indices for target onset, reach onset, touch screen, return and mouth
    targi = 1;
    TargetTimeix = COTData.Target{targi}.TargetTime;
    ReachTimeix = COTData.Target{targi}.ReachTimeix;
    TouchTimeix = COTData.Target{targi}.TouchTimeix;
    ReturnTimeix = COTData.Target{targi}.ReturnTimeix; 
    MouthTimeix = COTData.Target{targi}.MouthTimeix;
    [m, n] = size(TargetTimeix);
    if m == 1 || n == 1
        TargetTimeix = reshape(TargetTimeix, [m * n, 1]);
    end
    clear m n
    
    tbl_maTimeix = [table(TargetTimeix) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];
    
    % the tag of good reach trials
    mask_goodreach = COTData.Target{targi}.goodix_reach;
    % the tag of good return trials
    mask_goodreturn = COTData.Target{targi}.goodix_return;
    clear TargetTimeix ReachTimeix TouchTimeix ReturnTimeix MouthTimeix
    
    
    % varNames and varTypes for all return T_idxevent_ma, T_idxevent_lfp
    varNames_table = {'TargetTimeix', 'ReachTimeix', 'TouchTimeix', 'ReturnTimeix', 'MouthTimeix'};
    varTypes_table = {'double','double','double', 'double', 'double'};
    tbl_maTimeix.Properties.VariableNames = varNames_table;
    
    % total n_trial and n_event
    [n_trial, n_events] = size(tbl_maTimeix);
    
    % t_bef: time before target on, t_aft: time after mouth
    t_bef = 1; 
    t_aft = 0.5; 
    
    % extract T_idxevent_ma, smoothWspeed_trial, Wpos_smooth_trial, Wrist_smooth_trial
    n_bef = round(t_bef * fs_ma); % n_bef: index MA number before target on
    n_aft = round(t_aft * fs_ma); % n_aft: index MA number after mouth
    
    idx_strs = tbl_maTimeix.TargetTimeix - n_bef; % start ma index for each trial n_trial * 1
    idx_ends = tbl_maTimeix.MouthTimeix + n_aft; % end ma index for each trial n_trial * 1
    
    % smoothWspeed_trial: ntemp * ntrial
    smoothWspeed_trial = cell(n_trial, 1);
    Wpos_smooth_trial = cell(n_trial, 1);
    Wrist_smooth_trial = cell(n_trial, 1);
    
    for triali = 1:n_trial
        smoothWspeed_trial{triali} = COTData.smoothWspeed(idx_strs(triali) : idx_ends(triali));
        Wpos_smooth_trial{triali} = COTData.Wpos_smooth(idx_strs(triali) : idx_ends(triali));
        Wrist_smooth_trial{triali} = COTData.Wrist_smooth(idx_strs(triali) : idx_ends(triali), :);
    end
    T_idxevent_ma = table('Size', size(tbl_maTimeix), 'VariableTypes', varTypes_table, 'VariableNames',varNames_table);
    T_idxevent_ma{:, :} = tbl_maTimeix{:, :} - repmat(idx_strs, [1, n_events]);
    
    
    
    
    %% *_TDTwParadigm file for GoCue
    paradigmfilestruct = dir(fullfile(mafolder, '*COT_TDTwParadigm.mat'));
    if length(paradigmfilestruct) ~= 1
        
        disp([mafolder 'has ' num2str(length(paradigmfilestruct)) ' files, skip!'])
        
        lfptrial_cortical = []; lfptrial_dbs = [];
        mask_goodreach = []; mask_goodreturn = [];
        fs_lfp = []; T_idxevent_lfp = [];
        fs_ma = []; T_idxevent_ma = [];
        smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
        T_dbsChn = [];
        
        return
    end
    
    % load COTData from  *_TDTwParadigm.mat
    load(fullfile(mafolder, paradigmfilestruct.name), 'paradigm_data');
    StartTrial_TDT_ix = paradigm_data.TDTTrialInfo.StartTrial_TDT_ix;
    GoCue_TDT_ix = paradigm_data.TDTTrialInfo.GoCue_TDT_ix;
    ReachStart_TDT_ix = paradigm_data.TDTTrialInfo.ReachStart_TDT_ix;
    ReachEnd_TDT_ix = paradigm_data.TDTTrialInfo.ReachEnd_TDT_ix;
    EndTrial_TDT_ix = paradigm_data.TDTTrialInfo.EndTrial_TDT_ix;
    
    [m, n] = size(StartTrial_TDT_ix);
    if m == 1 || n == 1
        StartTrial_TDT_ix = reshape(StartTrial_TDT_ix, [m * n, 1]);
    end
    [m, n] = size(GoCue_TDT_ix);
    if m == 1 || n == 1
        GoCue_TDT_ix = reshape(GoCue_TDT_ix, [m * n, 1]);
    end
    [m, n] = size(ReachStart_TDT_ix);
    if m == 1 || n == 1
        ReachStart_TDT_ix = reshape(ReachStart_TDT_ix, [m * n, 1]);
    end
    [m, n] = size(ReachEnd_TDT_ix);
    if m == 1 || n == 1
        ReachEnd_TDT_ix = reshape(ReachEnd_TDT_ix, [m * n, 1]);
    end
    [m, n] = size(EndTrial_TDT_ix);
    if m == 1 || n == 1
        EndTrial_TDT_ix = reshape(EndTrial_TDT_ix, [m * n, 1]);
    end
    clear m n
    
    mask_targ0 = (paradigm_data.TDTTrialInfo.PositionIndex == 0);
    StartTrial_TDT_ix = StartTrial_TDT_ix(mask_targ0);
    GoCue_TDT_ix = GoCue_TDT_ix(mask_targ0);
    ReachStart_TDT_ix = ReachStart_TDT_ix(mask_targ0);
    ReachEnd_TDT_ix = ReachEnd_TDT_ix(mask_targ0);
    EndTrial_TDT_ix = EndTrial_TDT_ix(mask_targ0);
    clear mask_targ0

    
    %% LFP data
    % read each channel data in  LFP data to lfpdata (initialized when is the 1st channel)
    folder_cortical = fullfile(onedaypath, 'LFP', ['Block-' num2str(tdtbk)]);

    % files in folder_cortical are empty
    if isempty(dir(folder_cortical))
        disp([folder_cortical ' has no files!'])
        
        lfptrial_cortical = []; lfptrial_dbs = [];
        mask_goodreach = []; mask_goodreturn = [];
        fs_lfp = []; T_idxevent_lfp = [];
        fs_ma = []; T_idxevent_ma = [];
        smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
        T_dbsChn = [];
        
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
    for i = 1:length(chns)
        filename = [file_prefix num2str(chns(i)) '.nex'];
        [nexlfp_cortical] = readNexFile(fullfile(folder_cortical, filename));

        % extract the number of the structure containing LFPchn* data
        name_list = extractfield(cell2mat(nexlfp_cortical.contvars), 'name');
        i_lfp = find(contains(name_list, 'LFP')); %% i.e nexlfp_utah.contvars name == 'LFPch1', or 'MUAch1'

        if i == 1% first channel
            fs_lfpcortical = nexlfp_cortical.contvars{i_lfp}.ADFrequency;
            
            % add GoCueIndex for T_idxevent_ma, timeixtbl_ma 
            fs_TDT = nexlfp_cortical.freq; % used for GoCue_TDT_ix
            
            StartTrialix = round(StartTrial_TDT_ix / fs_TDT * fs_ma);
            GoCueix = round(GoCue_TDT_ix / fs_TDT * fs_ma);
            ReachStartix = round(ReachStart_TDT_ix / fs_TDT * fs_ma);
            ReachEndix = round(ReachEnd_TDT_ix / fs_TDT * fs_ma);
            EndTrialix = round(EndTrial_TDT_ix / fs_TDT * fs_ma);
            
            tbl_maFromTDTix = table(StartTrialix, GoCueix, ReachStartix, ReachEndix, EndTrialix);

            % the time index in the LFP neural data based on MA data
            tbl_lfpTimeix = tbl_maTimeix;
            tbl_lfpTimeix{:, :} = round(tbl_maTimeix{:, :} / fs_ma * fs_lfpcortical);

            % initialize lfp_utah : chn_lfp * n_temporal * n_trial
            n_bef = round(t_bef * fs_lfpcortical); % n_bef: index number before target on
            n_aft = round(t_aft * fs_lfpcortical); % n_aft: index number after mouth
            idx_str = tbl_lfpTimeix.TargetTimeix - n_bef;
            idx_end = tbl_lfpTimeix.MouthTimeix + n_aft;
            lfptrial_cortical = cell(n_trial,1);

            % idxtbl_lfptrialutah: the idx for events of target onset, reach onset, touch screen,
            % return and mouth in the trial matrix lfpdata_utah (chn_lfputah * n_temporal * n_trial)
            % the first sample corresponds to target onset - t_bef
            T_idxevent_lfp = tbl_lfpTimeix;
            T_idxevent_lfp{:, :} = T_idxevent_lfp{:, :} - repmat(idx_str, [1, n_events]);

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
                disp(['idx in ma' num2str(tbl_maTimeix{1, end})])
                disp(size(nexlfp_cortical.contvars{i_lfp}.data, 1))
                disp(idx_end(triali))
                disp(['maxlen  = ' num2str(maxlen)])
            end

            lfptrial_cortical{triali, 1} = cat(2, lfptrial_cortical{triali, 1}, nexlfp_cortical.contvars{i_lfp}.data(idx_str(triali):idx_end(triali)));

        end

        clear filename i_lfp
    end


    %% DBSLFP data
    % read DBSLFP in nex file
    dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(tdtbk)]); % 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\DBSLFP\Block-1'

    dbsfilepattern = fullfile(dbslfpfolder, '*DBSLFP.nex');
    dbsfiles = dir(dbsfilepattern);
    
    % dbs file does not exist
    if length(dbsfiles) ~= 1
        disp([dbsfilepattern ' has ' num2str(length(dbsfiles)) ' file, skip!'])
        
        lfptrial_cortical = []; lfptrial_dbs = [];
        mask_goodreach = []; mask_goodreturn = [];
        fs_lfp = []; T_idxevent_lfp = [];
        fs_ma = []; T_idxevent_ma = [];
        smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
        T_dbsChn = [];

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
        
        lfptrial_cortical = []; lfptrial_dbs = [];
        mask_goodreach = []; mask_goodreturn = [];
        fs_lfp = []; T_idxevent_lfp = [];
        fs_ma = []; T_idxevent_ma = [];
        smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
        T_dbsChn = [];
        
        return
    end

    if  ~(abs(convars(idx_dbs(1)).ADFrequency - fs_lfpcortical) < 0.001)
        % check the sampling frequency of dbs is the same as the lfp stored in separate channels or not
        disp(['the sampling frequency of dbs is not the same as the lfp stored in separate channels'])
        
        lfptrial_cortical = []; lfptrial_dbs = [];
        mask_goodreach = []; mask_goodreturn = [];
        fs_lfp = []; T_idxevent_lfp = [];
        fs_ma = []; T_idxevent_ma = [];
        smoothWspeed_trial = []; Wpos_smooth_trial = []; Wrist_smooth_trial = [];
        T_dbsChn = [];
        
        return
    end

    fs_lfp = fs_lfpcortical; % sampling frequency for dbs channels
    clear idxtbl_lfptrial_sepchn

    % extract each trial for lfp dbs data
    nchn_dbs = length(idx_dbs);
    lfptrial_dbs = cell(n_trial,1);

    for i = 1:nchn_dbs
        chni = idx_dbs(i);
        lfp_1chn = convars(chni).data;

        for triali = 1:n_trial
            lfptrial_dbs{triali, 1} = cat(2, lfptrial_dbs{triali, 1}, lfp_1chn(idx_str(triali):idx_end(triali)));
        end

        clear chni lfp_1chn triali
    end

    % chantbl_dbs: a table describes each dbs channel information
    elecchn = extractfield(convars(idx_dbs), 'name');
    elecchn = elecchn';
    area = cell(nchn_dbs, 1);
    area(1:8) = {'STN'};
    area(9:16) = {'GP'};
    T_dbsChn = [table(area) table(elecchn)];
    clear varName

    clear convars filename idx_stn idx_gp
end