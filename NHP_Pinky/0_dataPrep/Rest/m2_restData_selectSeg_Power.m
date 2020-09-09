function m2_restData_selectSeg_Power()
    %   Manually marked the good and bad segments
    %
    % 1. remove the segment manually marked bad ('n')
    %
    % 2. bipolar for STN and GP channels
    %
    % 3. Down sample trials into fs_new = 500
    %
    % 4. add variable GMChnAreas from GMChnsarea file
    %

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

    %% input setup
    inputfolder = fullfile(codecorresParentfolder, 'm1_restData_cleaned_extract');

    twin = 2;
    toverlap = twin * 0.9;
    freqs_roi = [5 40];
    fs_new = 500;

    animal = 'Pinky';

    % GrayMatter chn-area information file
    filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
    file_GMChnsarea = fullfile(datafolder, 'Pinky', filename_GMChnsarea);

    %% save setup
    savefolder = codecorresfolder;

    %% Start Here

    % GM chn Areas table
    nGM = 32;
    T_chnsarea_GM = chanInf_GM(file_GMChnsarea, nGM);

    % segment check
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);

    for fi = 1:nfiles
        file = fullfile(files(fi).folder, files(fi).name);

        load(file, 'fs', 'data_segments', 'chans_m1')

        disp([num2str(fi) '/' num2str(nfiles) ' ' files(fi).name]);

        % nwin and noverlap for pwelch calculation
        nwin = round(twin * fs);
        noverlap = round(toverlap * fs);

        nsegs = length(data_segments);
        segsBad = [];
        % pwelch for each segi
        for segi = 1:nsegs

            % lfp_m1: ntemp * nchns
            lfp_m1 = data_segments(segi).lfp_m1;
            lfp_meanm1 = mean(lfp_m1, 2);

            [pxx_m1, F_m1] = pwelch(lfp_meanm1, nwin, noverlap, [freqs_roi(1):1 / twin:freqs_roi(2)], fs);

            % lfp_stn: ntemp * 8; lfp_diffstn: ntemp * 7
            lfp_stn = data_segments(segi).lfp_stn;
            lfp_diffstn = diff(lfp_stn, [], 2);

            [pxx_stn, F_stn] = pwelch(lfp_diffstn, nwin, noverlap, [freqs_roi(1):1 / twin:freqs_roi(2)], fs);

            % lfp_gp: ntemp * 8; lfp_diffgp: ntemp * 7
            lfp_gp = data_segments(segi).lfp_gp;
            lfp_diffgp = diff(lfp_gp, [], 2);

            [pxx_gp, F_gp] = pwelch(lfp_diffgp, nwin, noverlap, [freqs_roi(1):1 / twin:freqs_roi(2)], fs);

            % plot pxx_m1, pxx_stn and pxx_gp
            figure('units', 'normalized', 'outerposition', [0.1 0. 0.75 0.75])

            subplot(4, 4, 1); plot(F_m1, pxx_m1); legend('m1'); xlim(freqs_roi);

            for chi = 1:size(lfp_diffstn, 2)
                subplot(4, 4, chi + 1); plot(F_stn, pxx_stn(:, chi)); legend(['stn - ' num2str(chi)]);
                xlim(freqs_roi);
            end

            for chi = 1:size(lfp_diffstn, 2)
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
                segsBad = [segsBad;segi];
            end

            %%% bipolar stn and gp %%%%
            lfp_stn = lfp_diffstn;
            lfp_gp = lfp_diffgp;

            %%% down sample %%%
            lfp_m1 = resample(lfp_m1, round(fs_new), round(fs));
            lfp_GM = resample(data_segments(segi).lfp_GM, round(fs_new), round(fs));
            lfp_stn = resample(lfp_stn, round(fs_new), round(fs));
            lfp_gp = resample(lfp_gp, round(fs_new), round(fs));
            

            %%% combine lfp_m1, lfp_stn and lfp_gp %%%%
            data_segments(segi).lfp = cat(2, lfp_m1, lfp_GM, lfp_stn, lfp_gp);

            %%% add chn-area information T_chnsarea  %%%
            if ~exist('T_chnsarea', 'var')
                nM1 = size(lfp_m1, 2);
                nGM = size(lfp_GM, 2);
                nSTN = size(lfp_stn, 2);
                nGP = size(lfp_gp, 2);

                T_chnsarea_M1 = chanInf_M1(chans_m1);
                T_chnsarea_DBS = chanInf_DBS(nSTN, nGP);

                T_chnsarea_GM.chni = T_chnsarea_GM.chni + nM1; 
                T_chnsarea_DBS.chni =  T_chnsarea_DBS.chni + nM1 + nGM;
                T_chnsarea = [T_chnsarea_M1; T_chnsarea_GM; T_chnsarea_DBS];

                clear nM1 nGM nSTN nGP T_chnsarea_M1 T_chnsarea_DBS
            end

            close all
            clear lfp_m1 lfp_mean pxx reply

        end

        % new fs
        fs = fs_new;

        % remove the bad segments
        data_segments(segsBad) = [];

        % remove lfp_m1, lfp_GM, lfp_stn and lfp_gp as they have been combined together
        data_segments = rmfield(data_segments, {'lfp_m1', 'lfp_GM', 'lfp_stn', 'lfp_gp'});

        savefile = fullfile(savefolder, files(fi).name);
        save(savefile, 'fs', 'data_segments', 'T_chnsarea');

        clear file chans_m1 fs data_segments
        clear nwin noverlap segi nsegs segsBad
    end

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
    T_chnsarea.chni = uint8([1: nGM]');
    T_chnsarea.brainarea = GMChnAreas;
    T_chnsarea.recordingchn = uint8([1: nGM]') + chi_firstGM -1;
    T_chnsarea.electype = cell(nGM,1);
    T_chnsarea.electype(:) = {'Gray Matter'};
    T_chnsarea.notes = cell(nGM,1);

end

function T_chnsarea = chanInf_M1(chans_m1)
    % extract M1 channel inf table
    %   Args:
    %       chans_m1: a vector containing m1 channel numbers
    %
    %   Output:
    %       T_chnsarea: table of M1 channel inf,
    %                   (T_chnsarea.Properties.VariableNames:
    %                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

    nM1 = length(chans_m1);
    chans_m1 = reshape(chans_m1, nM1, 1);
    

    % channel information table of M1
    T_chnsarea = table;
    T_chnsarea.chni = uint8([1: nM1]');
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

    T_chnsarea.chni = uint8([1: nSTN + nGP]');

    T_chnsarea.brainarea = cell(nSTN + nGP, 1);
    T_chnsarea.brainarea(1:nSTN) = {'STN'}; T_chnsarea.brainarea(nSTN + 1: nSTN + nGP) = {'GP'};

    T_chnsarea.recordingchn = [uint8([1: nSTN]');uint8([1: nGP]')];

    T_chnsarea.electype = cell(nSTN + nGP, 1);
    T_chnsarea.electype(:) = {'DBS'};

    T_chnsarea.notes = cell(nSTN + nGP, 1);

end