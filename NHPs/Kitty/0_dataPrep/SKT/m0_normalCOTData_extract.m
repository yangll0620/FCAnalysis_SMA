function m0_normalCOTData_extract()    
    %% extract the COT normal data for Kitty 
    %
    %	Inputs:
    %
    %		/home/lingling/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/data/Kitty/KittyNormal_SKB_ArrayDBSLFP_ALL.mat
    %
    %   Steps:
    %       1. extract trials from Data from Ying
    %           a. trial length = max(each trial length) + t_bef + t_aft 
    %                   t_bef: the time before target on (default: t_bef = 1)
    %                   t_aft: the time after mouth (default: t_aft = 0.5)
    %           b. only remain the trials markedin both good reach and good return
    %
    %
    %       2. bipolar for DBS channels (already done in input file)
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


    % datafolder, pipelinefolder
    [datafolder, ~, ~, ~] = exp_subfolders();
    correspipelinefolder = code_corresfolder(codefilepath, true, false);

    %% input setup
    tmpdata = regexp(codefilepath, 'NHPs/[a-zA-Z]*', 'match');
    animal = tmpdata{1}(length('NHPs/')+1:end);
    clear tmp
    
    % I think is COT other than SKB even if Ying named it SKB
    input_COTfile = fullfile(datafolder, animal, 'KittyNormal_SKB_ArrayDBSLFP_ALL.mat');
    
    
    % channel number of utah array, grayMatter, and dbs leads
    nM1 = 1; nSTN = 7; nGP = 7;
    
    % t_bef: time before target on, t_aft: time after mouth
    t_bef = 1; t_aft = 0.5; 
    
    % downsample
    fs_new = 500;
    
    pdcondition = 'normal';
    
   
    strformat_date_mat = 'mmddyy'; % the format of date string in mat, e.g. '101614'

    %% save setup
    savefolder = correspipelinefolder;
    savefilename_prefix = [animal '_COTData_'];

    strformat_date_save = 'yyyymmdd'; % the format of date string in saved filename, e.g. '20170123'
    
    %% Starting Here
    load(input_COTfile, 'data');
    
    %%% add chn-area information T_chnsarea  %%%
    T_chnsarea = chanInf_M1DBS(nM1, nSTN, nGP);
    
    %%% extract all trials for each day %%%
    nrecords = length(data);
    for ni = 1:nrecords
        tmpdata = data(ni);
        
        % extract goodind
        goodind = tmpdata.goodind{1};
        for i = 2: length(tmpdata.goodind)
            goodind = intersect(goodind, tmpdata.goodind{i});
        end
        
        fs = tmpdata.lfpfs;
        
        eventInds = round(tmpdata.Eventtime(goodind, :) * fs);
        

        
        %%% extract lfptrial and idxeventtbl
        maxlen = max(eventInds(:, 5) - eventInds(:, 1)) + 1;  % maximum length across all trials (unit: ind)
        n_bef = round(t_bef * fs); % n_bef: index number before target on
        n_aft = round(t_aft * fs); % n_aft: index number after mouth
        n_temporal = maxlen + n_bef + n_aft;
        
        % start and end index for each trial
        idxs_str = eventInds(:, 1) - n_bef;
        idxs_end = eventInds(:, 1) + (maxlen -1) + n_aft;
        
        [n_trials, n_timevars] = size(eventInds);
        lfptrial_M1 = zeros(nM1, n_temporal, n_trials);
        lfptrial_DBS = zeros(nSTN + nGP, n_temporal, n_trials);
        for triali = 1: size(eventInds, 1)
            lfptrial_M1(:, :, triali) = transpose(tmpdata.Arraylfp_mean(idxs_str(triali):idxs_end(triali), :));
            lfptrial_DBS(:, :, triali) = transpose(tmpdata.DBSlfp(idxs_str(triali):idxs_end(triali), :));
        end
        
        % cat M1 and DBS LFP %%%
        lfptrial = cat(1, lfptrial_M1, lfptrial_DBS);
          
        % the time index
        varNames =  {'TargetTimeix','ReachTimeix','TouchTimeix', 'ReturnTimeix', 'MouthTimeix'};
        T_idxevent = table;
        for evi = 1 : n_timevars
            T_idxevent = [T_idxevent table(eventInds(:, evi) - idxs_str, 'VariableNames', varNames(evi))];
        end
        
        
        
        
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
        dateofexp = datenum(tmpdata.Foldername(length(animal)+2:end), strformat_date_mat);
        
        
        % save
        savefilename = [savefilename_prefix pdcondition '_' datestr(dateofexp, strformat_date_save) '_n' num2str(ni) '.mat'];
        save(fullfile(savefolder, savefilename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');


        
        
        clear tmpdata
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




