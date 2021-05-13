function m2_restData_selectSeg_Power()
    %   Manually marked the good and bad segments
    %
    %
    %   Processing steps as follows:
    %       1. remove the segment manually marked bad ('n')
    %
    %       2. bipolar for STN and GP channels
    %
    %       3. Down sample trials into fs_new = 500
    %
    %       4. combine lfp_GM, lfp_stn and lfp_gp into lfp(ntemp * nchns)
    %
    %       5. add table T_chnsarea (nchns * 1  cell)
    %               For Bug: depth of channels was added
    %               if the depth record for one date  doesn't exist, skip that day

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
    
    
    [datafolder, ~, ~, ~] = exp_subfolders();
    
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
    inputfolder = fullfile(codecorresParentfolder, 'm1_restData_cleaned_extract');
    
    file_GMChnsarea = fullfile(datafolder, animal, [animal '_GMChnAreaInf.csv']);
    file_chnDepth = fullfile(codecorresParentfolder, '..', 'm1_dailyDepth', 'Bug_dailyDepth.xlsx');
    
    twin = 2;
    toverlap = twin * 0.9;
    freqs_roi = [5 40];

    refchns_m1 = [56 58 79 69 65];

    fs_new = 500;

    %% save setup
    savefolder = codecorresfolder;

    %% Start Here

    % segment check
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);

    for fi = 1:nfiles
        
        filename = files(fi).name;
        file = fullfile(files(fi).folder, files(fi).name);
        
        [i, j] = regexp(filename, '[0-9]*_tdt');
        dateofexp = datenum(filename(i:j- length('_tdt')), 'yyyymmdd');
        clear i j

        load(file, 'fs', 'data_segments')
        
        
        %%% add chn-area information T_chnsarea  %%%
        nGM = size(data_segments(1).lfp_array, 2);
        nSTN = size(data_segments(1).lfp_stn, 2);
        nGP = size(data_segments(1).lfp_stn, 2);

        if contains(filename, '_normal_')
            pdcond = 'normal';
        end
        if contains(filename, '_mild_')
            pdcond = 'mild';
        end
        if contains(filename, '_moderate_')
            pdcond = 'moderate';
        end

        T_chnsarea_GM = chanInf_GM(file_GMChnsarea, nGM, pdcond);
        T_chnsarea_GM_Depth = add_daily_GMDepth(T_chnsarea_GM, file_chnDepth, dateofexp, animal);
        if isempty(T_chnsarea_GM_Depth)
            clear nGM nSTN nGP  T_chnsarea_GM T_chnsarea_GM_Depth
            continue;
        end

        T_chnsarea_DBS = chanInf_DBS(nSTN - 1, nGP - 1);
        T_chnsarea_DBS_Depth = addvars(T_chnsarea_DBS, NaN(height(T_chnsarea_DBS), 1), 'NewVariableNames', 'depth', 'After', 'recordingchn');

        T_chnsarea = [T_chnsarea_GM_Depth; T_chnsarea_DBS_Depth];
        T_chnsarea.chni = [1: height(T_chnsarea)]';

        clear nGM nSTN nGP pdcondition T_chnsarea_GM T_chnsarea_GM_Depth T_chnsarea_DBS T_chnsarea_DBS_Depth


        disp([num2str(fi) '/' num2str(nfiles) ' ' filename]);

        % nwin and noverlap for pwelch calculation
        nwin = round(twin * fs);
        noverlap = round(toverlap * fs);

        %%% manually select segment %%%%
        nsegs = length(data_segments);
        segsBad = [];
        % pwelch for each segi
        for segi = 1:nsegs

            % lfp_array: ntemp * nchns
            lfp_array = data_segments(segi).lfp_array;
            lfp_meanm1 = mean(lfp_array(:, refchns_m1), 2);

            [pxx_m1, F_m1] = pwelch(lfp_meanm1, nwin, noverlap, [freqs_roi(1):1 / twin:freqs_roi(2)], fs);

            % lfp_stn: ntemp * 8; lfp_diffstn: ntemp * 7
            lfp_stn = data_segments(segi).lfp_stn;
            lfp_diffstn = diff(lfp_stn, [], 2);

            [pxx_stn, F_stn] = pwelch(lfp_diffstn, nwin, noverlap, [freqs_roi(1):1 / twin:freqs_roi(2)], fs);

            % lfp_gp: ntemp * 8; lfp_diffgp: ntemp * 7
            lfp_gp = data_segments(segi).lfp_gp;
            lfp_diffgp = diff(lfp_gp, [], 2);

            [pxx_gp, F_gp] = pwelch(lfp_diffgp, nwin, noverlap, [freqs_roi(1):1 / twin:freqs_roi(2)], fs);

            % figure
            figure('units', 'normalized', 'outerposition', [0.1 0. 0.75 0.75])
            subplot(4, 4, 1); plot(F_m1, pxx_m1); legend('m1'); xlim(freqs_roi);

            for chi = 1:size(lfp_diffstn, 2)
                subplot(4, 4, chi + 1); plot(F_stn, pxx_stn(:, chi)); legend(['stn - ' num2str(chi)]);
                xlim(freqs_roi);
            end

            for chi = 1:size(lfp_diffgp, 2)
                subplot(4, 4, chi + 8); plot(F_gp, pxx_gp(:, chi)); legend(['gp - ' num2str(chi)]);
                xlim(freqs_roi);
            end

            % manualy check the good('y') and the bad (n) segment
            reply = input(['segi=' num2str(segi) '/' num2str(nsegs) ', remain this seg? y/n [y]:'], 's');

            if isempty(reply)
                reply = 'y';
            end

            reply = lower(reply);

            if reply == 'n'
                segsBad = [segsBad; segi];
                close all
                continue;
            end

            %%% bipolar stn and gp %%%%
            lfp_stn = lfp_diffstn;
            lfp_gp = lfp_diffgp;

            %%% down sample %%%
            lfp_array = resample(lfp_array, round(fs_new), round(fs));
            lfp_stn = resample(lfp_stn, round(fs_new), round(fs));
            lfp_gp = resample(lfp_gp, round(fs_new), round(fs));

            %%% combine lfp_array, lfp_stn and lfp_gp %%%%
            data_segments(segi).lfp = cat(2, lfp_array, lfp_stn, lfp_gp);

            close all
            
            clear lfp_array lfp_meanm1 pxx_m1 F_m1
            clear lfp_stn lfp_diffstn pxx_stn F_stn lfp_gp lfp_diffgp pxx_gp F_gp
            clear reply

        end
        
        
        % new fs
        fs = fs_new;

        % remove the bad segments
        data_segments(segsBad) = [];

        % remove lfp_array, lfp_stn and lfp_gp as they have been combined together
        data_segments = rmfield(data_segments, {'lfp_array', 'lfp_stn', 'lfp_gp'});

        savefile = fullfile(savefolder, files(fi).name);
        save(savefile, 'fs', 'data_segments', 'T_chnsarea');

        clear file fs data_segments T_chnsarea savefile
        clear nwin noverlap segi nsegs segsBad
    end

end

function T_chnsarea = chanInf_GM(file_GMChnsarea, nGM, pdcondition)
    % extract gray matter channel inf table
    %
    %   Args:
    %       file_GMCchnsarea: the file storing the gray matter chn-area inf
    %       nGM: the total number of gray matter channels
    %       pdcondition: pd condition
    %
    %   Output:
    %       T_chnsarea: table of Gray matter channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    % initial the attitudes of chanarea table T_chnsarea:
    % chni_vec, electype, brainareas, notes, recordingchn
    chni_vec = uint8([1:nGM]');
    electype = cell(nGM, 1);
    brainareas = cell(nGM, 1);
    notes = cell(nGM, 1);
    recordingchn = zeros(nGM, 1);
    electype(1:nGM, 1) = {'Gray Matter'}; % electype

    % deal with Gray Matter
    brainareas(1:nGM, 1) = {''};
    T = readtable(file_GMChnsarea);

    if strcmp(pdcondition, 'normal')
        T.channels = T.channels_normal;
    end

    if strcmp(pdcondition, 'mild')
        T.channels = T.channels_mild;
    end

    if strcmp(pdcondition, 'moderate')
        T.channels = T.channels_moderate;
    end

    for i = 1:length(T.brainarea)
        area = T.brainarea{i};
        tmpcell = split(T.channels{i}, ',');

        for j = 1:length(tmpcell)
            chn = str2num(char(tmpcell{j}));
            brainareas(chn, 1) = {area};

            clear chn
        end

    end

    recordingchn(1:nGM) = [1:nGM]';
    notes(1:nGM, 1) = {''};

    % channel information table
    T_chnsarea = table;
    T_chnsarea.chni = chni_vec;
    T_chnsarea.brainarea = brainareas;
    T_chnsarea.recordingchn = recordingchn;
    T_chnsarea.electype = electype;
    T_chnsarea.notes = notes;
    
    
    % remove the brainarea: VA/Vlo/STN, Gpe/i
    mask_remove = cellfun(@(x) contains(x, '/'), T_chnsarea.brainarea);
    T_chnsarea{mask_remove, 'brainarea'} = {''};

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

    T_chnsarea.chni = uint8([1: nSTN + nGP]');

    T_chnsarea.brainarea = cell(nSTN + nGP, 1);
    T_chnsarea.brainarea(1:nSTN) = {'STN'}; T_chnsarea.brainarea(nSTN + 1: nSTN + nGP) = {'GP'};

    T_chnsarea.recordingchn = [uint8([1: nSTN]');uint8([1: nGP]')];

    T_chnsarea.electype = cell(nSTN + nGP, 1);
    T_chnsarea.electype(:) = {'DBS'};

    T_chnsarea.notes = cell(nSTN + nGP, 1);

end


function T_chnsarea_GM_Depth = add_daily_GMDepth(T_chnsarea_GM, file_chnDepth, dateofexp, animal)
    %%
    %   Args:
    %       T_chnsarea_GM: table for GM channels area without depth variable
    %       file_chnDepth: daily depth file containing normal, mild and moderate sheets(e.g. Bug_dailyDepth.xlsx)
    %       dateofexp: a particular exp date in datenum format (e.g datenum('04022019', 'mmddyyyy'))
    %
    %   return:
    %       T_chnsarea_GM_Depth: table for GM channels area with depth variable
    %                            return [] if the depthh record on dateofexp doesn't exist
    %                            the depth of chni will be NAN if brainarea = ''

    %% code Start here
    
    if nargin < 4
        animal = 'Bug';
    end
    strformat_date = 'dd-mmm-yyyy'; % the format of date string in folder in root2, e.g '012317'
    pdcond = parsePDCondition(dateofexp, animal);
    
    
    tbl = readtable(file_chnDepth, 'sheet', pdcond);

    % separate into t_chnDepth and t_chnArea
    t_chnDepth = tbl(2:end, :);
    t_chnArea = tbl(1, :);



    %%% replace date with a datenum variable and sort t_depth based on datenum %%%
    datenums = cell2mat(cellfun(@(x) datenum(x, strformat_date), t_chnDepth.Date,'UniformOutput', false));
    t_chnDepth = addvars(t_chnDepth, datenums, 'After', 'Date', 'NewVariableNames', 'datenum');
    t_chnDepth = removevars(t_chnDepth, 'Date');
    t_chnDepth = sortrows(t_chnDepth, {'datenum'});



    %%% add and extract new variable depth for T_chnsarea_GM  %%%

    % find the idx_dateDepth for dateofexp
    idx_dateDepth = find(t_chnDepth.datenum == dateofexp);
    if isempty(idx_dateDepth)  % if no equal, find the date before
        disp(['no depth record for' datestr(dateofexp)])
        T_chnsarea_GM_Depth = [];
        return
    end

    if length(idx_dateDepth) > 1
        disp(['length of dateDepth = ' num2str(length(idx_dateDepth))]);
        
        T_chnsarea_GM_Depth = [];
        return;
    end


    % assign each  channel with its depth on dateofexp if has
    T_chnsarea_GM_Depth = addvars(T_chnsarea_GM, NaN(height(T_chnsarea_GM), 1), 'NewVariableNames', 'depth', 'After', 'recordingchn');
    for vari = 2 : width( t_chnDepth)

        varName = t_chnDepth.Properties.VariableNames{vari}; % varNames: chan4, chan97 
        chn = str2num(varName(length('chan')+1:end));
        
        idx = find(T_chnsarea_GM_Depth.recordingchn == chn);
        
        if isempty(T_chnsarea_GM_Depth.brainarea{idx}) % empty brain area in T_chnsarea_GM_Depth
            continue;
        end

        if ~strcmp(t_chnArea{1, vari}, T_chnsarea_GM_Depth.brainarea{idx})
            disp(['brain area not equal for ' num2str(chn)])
            continue;
        end


        T_chnsarea_GM_Depth{idx, 'depth'} = str2num(cell2mat(t_chnDepth{idx_dateDepth, vari}));

        clear tmp area1 area2
        clear varName chn idx
    end
end