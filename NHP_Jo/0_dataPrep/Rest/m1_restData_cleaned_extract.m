function m1_restData_cleaned_extract()
    %% extract cleaned restData using files Ying used
    %
    %   1. remove data of channels in m1array_to_remove
    %
    %   2. remove data with eye-closed
    %
    %   3. remove data marked with movement
    %
    %   Inputs:
    %       configFile: datafolder/ config_m1lf_fromYing.mat
    %
    %       dateBlocksXLSFile: ../m0_restData_dateBlockYingUsed/dateBlocksYingUsed_rest.xlsx
    %
    %       filesUsed:  m0_restData_uniqdatesYingUsed_extract/uniqFileYingUsed_Pinky.mat


    %% folder generate
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');

    % find the codefolder
    idx = strfind(codefilepath, 'code');
    codefolder = codefilepath(1:idx + length('code') - 1);
    clear idx

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));

    % the corresponding pipeline folder for this code
    [codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

    % datafolder
    [datafolder, ~, ~, ~] = exp_subfolders();

    %% global parameter
    animal = 'Jo';

    % to be removed channels in m1, from Ying
    m1array_to_remove = [10 7; 10 8; 10 9; 10 10];

    %% input setup

    % input folder: root2 in server
    folder_processed_root2 = fullfile('/home/lingling/root2/Animals', animal, 'Recording/Processed/DataDatabase');

    % threshold used by Ying to extract cleaned resting data
    configFile = fullfile(datafolder, 'config_m1lf_fromYing.mat');

    % dateBlock File
    dateBlocksXLSFile = fullfile(codecorresParentfolder, 'm0_restData_dateBlockYingUsed', 'dateBlocksYingUsed_rest.xlsx');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_prefix = [animal '_cleanedRestData_'];
    savename_datestr_format = 'yyyymmdd';

    %% global variables

    % start time point, data before t0 are removed
    t0 = 3;

    % only remain the interval whose duration is larger than min_comp_time
    min_comp_time = 15;

    % threshold for eye_open
    thr_eye_open = 0.5000;

    %%  Starting here

    % load config_m1lf, and dateBlocks
    load(configFile, 'config_m1lf');
    tbl_dateBlocks = readtable(dateBlocksXLSFile);

    nfiles = height(tbl_dateBlocks);

    for i = 1:nfiles
        dateblockstr = tbl_dateBlocks(i, :).dateBlock_rest{1}; % dateblockstr = "20151002_1"
        cond = tbl_dateBlocks(i, :).cond{1};

        % extract the date_num, and blocki
        tmp = split(dateblockstr, '_');
        date_num = datenum(tmp{1}, 'yyyymmdd');
        tdtblocki = str2num(tmp{2});

        if ~(strcmp(cond, 'normal') || strcmp(cond, 'mild'))
            continue;
        end

        folder_allsync = fullfile(folder_processed_root2, [animal '_' datestr(date_num, 'mmddyy')], ['Block-' num2str(tdtblocki)]);
        filepattern_allsync = fullfile(folder_allsync, [animal '*' datestr(date_num, 'yyyymmdd') '_*_all_sync.mat']);

        files_allsync = dir(filepattern_allsync);

        % file not exist
        if (length(files_allsync) ~= 1)
            disp([num2str(i) '/' num2str(nfiles) ':' filepattern_allsync ' has ' num2str(length(files_allsync)) ' files, not 1 file'])
            continue;
        end

        disp(['dealing ' num2str(i) '/' num2str(nfiles) ': ' fullfile(files_allsync.folder, files_allsync.name)])

        % load data
        load(fullfile(files_allsync.folder, files_allsync.name), 'data');

        % fs
        fs = data.lfp_stn_fs(1);

        if (fs < 3000)
            % fs < 3000, the designed filter in filtered_lfp() is not right
            continue;
        end

        % get the thr_power1
        if strcmp(data.condition, 'PD')
            eval(['thr_power1 = config_m1lf.max_' lower(animal) '_pd;'])
        else
            eval(['thr_power1 = config_m1lf.max_' lower(animal) '_normal;'])
        end

        [data_segments, segsIndex] = cleanedRestData(data, m1array_to_remove, thr_power1, thr_eye_open, t0, min_comp_time, config_m1lf.freq_band);

        if (isempty(data_segments))
            continue;
        end

        savefile = fullfile(savefolder, [savefilename_prefix cond '_' datestr(date_num, savename_datestr_format) '_tdt' num2str(tdtblocki) '.mat']);
        save(savefile, 'data_segments', 'segsIndex', 'fs')

        clear allsyncfile dateblockstr thr_power1
        clear tmp cond
    end

    disp('Done!')
end

function [data_segments, segsIndex] = cleanedRestData(data, m1array_to_remove, thr_power1, thr_eye_open, t0, min_comp_time, freq_band)

    %%
    fs = data.lfp_stn_fs(1);
    n0 = round(t0 * fs);
    min_samples = round(min_comp_time * fs);

    % get the lfp_m1 after removing the m1array_to_remove channels
    lfp_m1 = get_lfp_removeBadChns(data.lfp_array, m1array_to_remove);

    % get lfp_stn, lfp_gp, and lfp_GM data
    lfp_stn = zeros(length(data.lfp_stn{1}), length(data.lfp_stn));

    for chi = 1:length(data.lfp_stn)
        lfp_stn(:, chi) = data.lfp_stn{chi};
    end

    lfp_gp = zeros(length(data.lfp_gp{1}), length(data.lfp_gp));

    for chi = 1:length(data.lfp_gp)
        lfp_gp(:, chi) = data.lfp_gp{chi};
    end

    %%  get states for each time point of lfp, 1: remain this time point, 0: remove
    m1power_state_mintime = extract_m1power_state_mintime(mean(lfp_m1, 2), thr_power1, min_samples, fs, freq_band);

    ntemp = length(lfp_m1);
    time_us = [0:ntemp - 1] / fs;
    eye_state_mintime = extract_eye_state_mintime(data.eye_area_filt, thr_eye_open, time_us, min_samples);

    % awake state should be both eye is open and m1_power is less than thr_power1
    awake_state_mintime = ((eye_state_mintime >= 1) & (m1power_state_mintime >= 1));

    % start and end index pair for each no movement segments idx_NoMoveSegs: n_NoMoveSegs * 2 (idx_start, idx_end)
    idx_noMoveSegs = round(data.SponMovement.not_movement_times * fs);
    to_remove_ma = extract_to_remove_ma(idx_noMoveSegs, ntemp);

    states = awake_state_mintime & (~to_remove_ma);

    %% filter lfp_m1, lfp_stn, lfp_gp, and lfp_GM and truncate at n0
    [lfp_m1, lfp_stn, lfp_gp] = filtered_lfp(lfp_m1, lfp_stn, lfp_gp, fs);

    %% truncate
    lfp_m1(1:n0 - 1, :) = [];
    lfp_stn(1:n0 - 1, :) = [];
    lfp_gp(1:n0 - 1, :) = [];
    states(1:n0 - 1) = [];

    %% extract segments
    %  resting data segments indices vec_segsIndex: nsegs * 2
    segsIndex = get_segIndex(states, min_samples);

    % extract segments
    for seg_ind = 1:size(segsIndex, 1)

        idx_segStr = segsIndex(seg_ind, 1);
        idx_segEnd = segsIndex(seg_ind, 2);

        %  lfp_m1
        if ~exist('data_segments', 'var')
            data_segments(1).lfp_m1 = lfp_m1(idx_segStr:idx_segEnd, :);
        else
            data_segments(end + 1).lfp_m1 = lfp_m1(idx_segStr:idx_segEnd, :);
        end

        % lfp_stn and lfp_gp
        data_segments(end).lfp_stn = lfp_stn(idx_segStr:idx_segEnd, :);
        data_segments(end).lfp_gp = lfp_gp(idx_segStr:idx_segEnd, :);

        clear idx_segStr idx_segEnd
    end

    if ~exist('data_segments', 'var')
        data_segments = [];
        segsIndex = [];
    end

end

function [lfp_m1, lfp_stn, lfp_gp] = filtered_lfp(lfp_m1, lfp_stn, lfp_gp, fs)
    %% pass filter lfp_m1, lfp_stn, and lfp_gp

    % pass filters of m1, dbs and GM data
    [bhp, ahp] = butter(2, 2 * 2 / fs, 'high'); % high pass
    [blphf, alphf] = butter(2, 2 * 700 / fs); % low pass

    % high pass filter lfp_m1
    for chi = 1:size(lfp_m1, 2)
        lfp_m1(:, chi) = filtfilt(bhp, ahp, lfp_m1(:, chi));
    end

    % high pass filter lfp_stn and lfp_gp
    for chi = 1:size(lfp_stn, 2)
        % filter
        lfp_stn(:, chi) = filtfilt(bhp, ahp, lfp_stn(:, chi));
        lfp_gp(:, chi) = filtfilt(bhp, ahp, lfp_gp(:, chi));
    end

end

function to_remove_ma = extract_to_remove_ma(idx_noMoveSegs, n)
    %% extract the to_remove_ma with MA data

    noMovData = zeros(1, n);

    for i = 1:size(idx_noMoveSegs, 1)
        noMovData(idx_noMoveSegs(i, 1):idx_noMoveSegs(i, 2)) = 1;
    end

    to_remove_ma = ~noMovData;

end

function m1power_state_mintime = extract_m1power_state_mintime(lfp_m1, thr_power1, min_samples, fs, freq_band)
    %% extract m1power_state_mintime

    % extracting lfp_m1_env through band pass, abs(hilbert) and low pass
    lfp_m1_env = lfp_m1;

    % band pass
    [blf1, alf1] = butter(2, 2 * freq_band / fs);
    lfp_m1_env = filtfilt(blf1, alf1, lfp_m1_env);

    % hilbert
    lfp_m1_env = abs (hilbert(lfp_m1_env));

    % low pass f
    [blf2, alf2] = butter(2, 2 * 0.15 / fs);
    lfp_m1_env = filtfilt(blf2, alf2, lfp_m1_env);

    % set m1_power_state =0 for lfp_m1_env > thr_power1 , 1 for lfp_m1_env <= thr_power1
    m1_power_state = ones(1, length(lfp_m1_env)) * 0.5;
    m1_power_state(lfp_m1_env > thr_power1) = 0;
    m1_power_state(lfp_m1_env <= thr_power1) = 1;

    [m1power_state_mintime] = get_state_mintime(m1_power_state, min_samples);

end

function eye_state_mintime = extract_eye_state_mintime(eye_area_ds, thr_eye_open, time_us, min_samples)
    %% extract eye_state_mintime for each time point in based on eye_area_ds and thr_eye_open

    % time point for eye and lfp seperately
    time_ds = [0:length(eye_area_ds) - 1] / 15;

    % get eye_area with the same time resolution of lfp data using 1-D interpolation
    eye_area = interp1(time_ds, eye_area_ds, time_us);

    % set eye_state =0 for eye_area < thr_eye_open, 1 for eye_area >= thr_eye_open
    eye_state = ones(1, length(eye_area));
    eye_state(eye_area < thr_eye_open) = 0;
    eye_state(eye_area >= thr_eye_open) = 1;

    [eye_state_mintime] = get_state_mintime(eye_state, min_samples);

end
