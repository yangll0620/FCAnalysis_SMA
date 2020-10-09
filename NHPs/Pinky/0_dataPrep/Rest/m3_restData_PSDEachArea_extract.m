function m3_restData_PSDEachArea_extract()
    %%   PSD estimates for mild and normal individually
    %
    %       psd for each brain area, as well as each DBS contact
    %
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

    %%  input setup
    
    animal = 'Pinky';

    % pwelch psd estimate variable
    twin_pwelch = 2;

    % variables for plotting
    plotF_AOI = [5 50];

    % input folder
    inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_Power');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_prefix = 'psd_';

    %% Code Start Here

    file_psdall = fullfile(savefolder, [savefilename_prefix '_allsegs_normalmildmoderate.mat']);

    %%%  calculate/load dbs psd for mild and normal %%%
    if ~exist(file_psdall, 'file')% not exist

        [pxxs_allfiles_normal, F_pxx_normal] = pxx_eacharea_allfiles(dir(fullfile(inputfolder, '*_normal_*.mat')), twin_pwelch);
        [pxxs_allfiles_mild, F_pxx_mild] = pxx_eacharea_allfiles(dir(fullfile(inputfolder, '*_mild_*.mat')), twin_pwelch);
        [pxxs_allfiles_moderate, F_pxx_moderate] = pxx_eacharea_allfiles(dir(fullfile(inputfolder, '*_moderate_*.mat')), twin_pwelch);

        if ~isequal(F_pxx_normal, F_pxx_mild, F_pxx_moderate)
            disp('F_pxx_normal, F_pxx_mild and F_pxx_moderate not equal');

            return;
        else
            F_pxx = F_pxx_normal;

            save(file_psdall, 'pxxs_allfiles_normal', 'pxxs_allfiles_mild', 'pxxs_allfiles_moderate','F_pxx');
            clear F_pxx_normal F_pxx_mild
        end

    else
        load(file_psdall, 'pxxs_allfiles_normal', 'pxxs_allfiles_mild', 'pxxs_allfiles_moderate','F_pxx');
    end

    %%%  plot  %%%
    brainareas = fieldnames(pxxs_allfiles_mild);

    for i = 1:length(brainareas)
        brainarea = brainareas{i};

        % load normal and mild data
        eval(['psd_normal = pxxs_allfiles_normal.' brainarea ';'])
        eval(['psd_mild = pxxs_allfiles_mild.' brainarea ';'])
        eval(['psd_moderate = pxxs_allfiles_moderate.' brainarea ';'])

        if ~strcmp(brainarea, 'STN') &&~strcmp(brainarea, 'GP')
            % if not DBS case, psd_normal, psd_mild: nfs * nsegs
            plotPSD_comp_1chn(psd_normal, psd_mild, psd_moderate, F_pxx, plotF_AOI, savefolder, brainarea)

        else
            % if DBS case, pxxs_allfiles.GP, pxxs.GP: nfs * nchns * nsegs
            plotPSD_comp_multichns(psd_normal, psd_mild, psd_moderate, F_pxx, plotF_AOI, savefolder, brainarea)
        end

    end
    
    
    
    %%% combine all figures into one %%%%
    
    % empty figure
    figure
    set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
    annotation(gcf,'textbox',...
        [0.3 0.7 0.35 0.15],...
        'String',{[animal ' Rest Data']},...
        'LineStyle','none',...
        'FontSize',15, 'FontWeight','bold',...
        'FitBoxToText','off');
    saveas(gcf, fullfile(savefolder, 'text'), 'png')
    img_text =  imread(fullfile(savefolder, 'text.png'));
    
    
    %%% combine all figures into one %%%%
    close all
    
    %----DBS in one figure ---%
    brainarea = 'STN';
    imgs_col1 = []; imgs_col2 = []; % two columns
    for i = 1: 7
        img = imread(fullfile(savefolder,['psd_' brainarea '_ch' num2str(i) '.png']));
        
        if mod(i, 2) == 1
            imgs_col1 = cat(1, imgs_col1, img);
        else
            imgs_col2 = cat(1, imgs_col2, img);
        end
        
        clear img
    end
    imgs_col2 = cat(1, imgs_col2, img_text); % the last one in column2 is empty
    imgs_STN = cat(2, imgs_col1, imgs_col2);
    clear imgs_col1 imgs_col2
    
    brainarea = 'GP';
    imgs_col1 = []; imgs_col2 = []; % two columns
    for i = 1: 7
        img = imread(fullfile(savefolder,['psd_' brainarea '_ch' num2str(i) '.png']));
        
        if mod(i, 2) == 1
            imgs_col1 = cat(1, imgs_col1, img);
        else
            imgs_col2 = cat(1, imgs_col2, img);
        end
        
        clear img
    end
    imgs_col2 = cat(1, imgs_col2, zeros(size(img_text)) + 255); % the last one in column2 is empty
    imgs_GP= cat(2, imgs_col1, imgs_col2);
    
    
    imgs_DBS = cat(2, imgs_STN, imgs_GP);
    imwrite(imgs_DBS,  fullfile(savefolder, 'combinedDBS.png'));
    
    
    %----M1, SMA and Thalams in one figure ---%
    imgs_tha = [];
    thalamus = {'VA', 'VLo', 'VPLo'};
    for i = 1 : length(thalamus)
        brainarea = thalamus{i};

        imgl = imread(fullfile(savefolder,['psd_l' brainarea '.png'])); 
        imgr = imread(fullfile(savefolder,['psd_r' brainarea '.png'])); 
        
        img = cat(2, imgl, imgr);
        imgs_tha = cat(1, imgs_tha, img);
        
        clear imgl imgr img brainarea
    end
    
    img_M1 = imread(fullfile(savefolder,'psd_M1.png')); 
    img_lSMA = imread(fullfile(savefolder,'psd_lSMA.png')); 
    img_rSMA = imread(fullfile(savefolder,'psd_rSMA.png')); 
    imgs_SMA = cat(2, img_lSMA, img_rSMA);
    imgs_2M1 = cat(2, img_text, img_M1);
    imgs_M1SMA = cat(1, imgs_SMA, imgs_2M1, zeros(size(imgs_SMA)) + 255);
    clear img_M1 img_lSMA img_rSMA imgs_2M1 imgs_SMA
    
    img_M1SMATha = cat(2, imgs_tha, imgs_M1SMA);
    imwrite(img_M1SMATha,  fullfile(savefolder, 'combinedM1SMATha.png'));
    
end

function [pxxs_allfiles, F_pxx] = pxx_eacharea_allfiles(files, twin_pwelch)
    %% extract psd from all the files, each psd for each dbs contact and one psd for one area (except dbs)
    %
    % Arg:
    %       files = dir(fullfile('.', '*_mild_*.mat'));
    %       twin_pwelch : twin for pwelch
    %
    % Outputs:
    %
    %       pxxs_allfiles: the PSD estimate of all the segments from the files
    %           e.g. pxxs =
    %                   struct with fields:
    %                      M1: [nfs �� nsegs double]
    %                     STN: [nfs �� nSTNchns �� nsegs double]
    %                      GP: [nfs �� nGPchns �� nsegs double]
    %
    %       F_pxx: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)

    nfiles = length(files);

    for filei = 1:nfiles

        file = fullfile(files(filei).folder, files(filei).name);

        [pxxs_1file, F_pxx_1file] = pxx_eacharea_onefile(file, twin_pwelch);

        if (isempty(pxxs_1file))
            continue;
        end
        
        brainareas = fieldnames(pxxs_1file);

        
        % combined pxxs from all the files
        if (~exist('pxxs_allfiles', 'var'))
            pxxs_allfiles = pxxs_1file;
        else

            for i = 1:length(brainareas)
                brainarea = brainareas{i};

                if ~strcmp(brainarea, 'STN') &&~strcmp(brainarea, 'GP')
                    % if not DBS case, pxxs_allfiles.M1, pxxs.M1: nfs * nsegs
                    eval(['pxxs_allfiles.' brainarea '= cat(2, pxxs_allfiles.' brainarea ', pxxs_1file.' brainarea ');'])
                else
                    % if DBS case, pxxs_allfiles.GP, pxxs.GP: nfs * nchns * nsegs
                    eval(['pxxs_allfiles.' brainarea '= cat(3, pxxs_allfiles.' brainarea ', pxxs_1file.' brainarea ');'])
                end

                clear onearea
            end

        end

        % extract the F_pxx
        if ~exist('F_pxx', 'var')
            F_pxx = F_pxx_1file;
        else

            if F_pxx ~= F_pxx_1file
                disp(['F_pxx ~=f for ' brainarea ', segi = ' num2str(segi) '']);

                F_pxx = [];
                pxxs_allfiles = [];

                return;
            end

        end

        clear file pxxs_1file F_pxx_1file brainareas
    end

end

function [pxxs, F_pxx] = pxx_eacharea_onefile(file, twin_pwelch)
    %% extract psd of all the segments from the files, each psd for each dbs contact and one psd for one area (except dbs)
    %
    % Outputs:
    %
    %       pxxs: the PSD estimate of all the segments from the files
    %           e.g. pxxs =
    %                   struct with fields:
    %                      M1: [nfs * nsegs double]
    %                     STN: [nfs * nSTNchns * nsegs double]
    %                      GP: [nfs * nGPchns * nsegs double]
    %
    %       F_pxx: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)

    % load data
    load(file, 'fs', 'data_segments', 'T_chnsarea');

    if (isempty(data_segments))% data_segments is empty
        F_pxx = [];
        pxxs = [];
        return;
    end

    % extract uniqBrainAreas
    mask_emptyarea = cellfun(@(x) isempty(x), T_chnsarea.brainarea);
    uniqBrainAreas = unique(T_chnsarea.brainarea(~mask_emptyarea));

    % psd pwelch paramers
    nwins = round(twin_pwelch * fs);
    noverlap = round(nwins * 0.9);

    pxxs = struct();

    for i = 1:length(uniqBrainAreas)
        brainarea = uniqBrainAreas{i};

        % extract the idx of brainarea
        mask_area = strcmp(T_chnsarea.brainarea, brainarea);

        for segi = 1:length(data_segments)

            % extract the lfp data of segi in brainarea, lfp_oneseg: ntemp * nchns
            lfp_oneseg = data_segments(segi).lfp(:, mask_area);

            if ~strcmp(brainarea, 'STN') &&~strcmp(brainarea, 'GP')
                % if not DBS case, use averaged lfp
                lfp_oneseg = mean(lfp_oneseg, 2);
            end

            %%% calcualte PSD throuch zscore and pwelch %%%

            % zscore of the lfp_oneseg along each column
            lfp_zscore = zscore(lfp_oneseg);

            % PSD is computed independently for each column, Pxx: nfs * nchns, f: nfs * 1
            [Pxx, f] = pwelch(lfp_zscore, nwins, noverlap, nwins, fs);

            if ~exist('F_pxx', 'var')
                F_pxx = f;
            else

                if F_pxx ~= f
                    disp(['F_pxx ~=f for ' brainarea ', segi = ' num2str(segi) '']);

                    F_pxx = [];
                    pxxs = [];

                    return;
                end

            end

            if (~isfield(pxxs, brainarea))% first time calculate pxx for brainarea
                eval(['pxxs.' brainarea ' = Pxx;'])
            end

            if ~strcmp(brainarea, 'STN') &&~strcmp(brainarea, 'GP')
                % if not DBS case, Pxx: nfs * 1; pxxs.M1: nfs * nsegs
                eval(['pxxs.' brainarea ' = cat(2, pxxs.' brainarea ', Pxx);'])
            else
                % if DBS case, Pxx: nfs * nchns; pxxs.STN: nfs * nchns * nsegs
                eval(['pxxs.' brainarea ' = cat(3, pxxs.' brainarea ', Pxx);'])
            end

            clear lfp_oneseg lfp_zscore Pxx f
        end

        clear brainarea mask_area segi
    end

end


function plotPSD_comp_1chn(psd_normal, psd_mild, psd_moderate, F_all, plotF_AOI, savefolder, brainarea)
    %%  plot the psd comparison of normal and mild
    %
    %   Inputs:
    %       psd_normal, psd_mild, psd_moderate: psd of all segments in normal, mild and moderate, nfs * nsegs
    %
    %       F_all: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
    %
    %       plotF_AOI: the plot frequence range (in hertz)  (1 * 2)
    %
    %       savefile_prefix: the savefile prefix including folder

    % find the idx for F_AOI
    idx_AOI = find(F_all >= plotF_AOI(1) & F_all <= plotF_AOI(2));

    % colors setup
    color_normal_range = [224, 255, 255] / 255;
    color_normal_mean = [0, 0, 255] / 255;
    color_mild_range = [255, 228, 225] / 255;
    color_mild_mean = [255, 00, 0] / 255;
    color_moderate_range = [238, 238, 238] / 255;
    color_moderate_mean = [0, 0, 0] / 255;

    % plot setup
    linewidth = 1.5;

    % F_AOI
    F_AOI = F_all(idx_AOI);
    % reshape F_AOI
    [m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m * n); clear m n

    % extract the  psd of AOI
    psd_normal_FAOI = psd_normal(idx_AOI, :);
    psd_mild_FAOI = psd_mild(idx_AOI, :);
    psd_moderate_FAOI = psd_moderate(idx_AOI, :);

    psd_normal_high = max(psd_normal_FAOI, [], 2);
    psd_normal_low = min(psd_normal_FAOI, [], 2);
    psd_normal_mean = mean(psd_normal_FAOI, 2);

    psd_mild_high = max(psd_mild_FAOI, [], 2);
    psd_mild_low = min(psd_mild_FAOI, [], 2);
    psd_mild_mean = mean(psd_mild_FAOI, 2);


    psd_moderate_high = max(psd_moderate_FAOI, [], 2);
    psd_moderate_low = min(psd_moderate_FAOI, [], 2);
    psd_moderate_mean = mean(psd_moderate_FAOI, 2);

    % reshape, psd_*_high/low/mean into 1 * nfs

    [m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m * n); clear m n
    [m, n] = size(psd_normal_low); psd_normal_low = reshape(psd_normal_low, 1, m * n); clear m n
    [m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m * n); clear m n

    [m, n] = size(psd_mild_high); psd_mild_high = reshape(psd_mild_high, 1, m * n); clear m n
    [m, n] = size(psd_mild_low); psd_mild_low = reshape(psd_mild_low, 1, m * n); clear m n
    [m, n] = size(psd_mild_mean); psd_mild_mean = reshape(psd_mild_mean, 1, m * n); clear m n


    [m, n] = size(psd_moderate_high); psd_moderate_high = reshape(psd_moderate_high, 1, m * n); clear m n
    [m, n] = size(psd_moderate_low); psd_moderate_low = reshape(psd_moderate_low, 1, m * n); clear m n
    [m, n] = size(psd_moderate_mean); psd_moderate_mean = reshape(psd_moderate_mean, 1, m * n); clear m n

    % plot range
    figure
    set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
    fill([F_AOI flip(F_AOI)], [psd_normal_high flip(psd_normal_low)], color_normal_range, 'LineStyle', 'none')
    hold all
    fill([F_AOI flip(F_AOI)], [psd_mild_high flip(psd_mild_low)], color_mild_range, 'LineStyle', 'none')
    fill([F_AOI flip(F_AOI)], [psd_moderate_high flip(psd_moderate_low)], color_moderate_range, 'LineStyle', 'none')

    % plot mean
    h1 = plot(F_AOI, psd_normal_mean, 'Color', color_normal_mean, 'LineWidth', linewidth);
    h2 = plot(F_AOI, psd_mild_mean, 'Color', color_mild_mean, 'LineWidth', linewidth);
    h3 = plot(F_AOI, psd_moderate_mean, 'Color', color_moderate_mean, 'LineWidth', linewidth);

    % find the frequency with maximum density
    [maxPSD, idx_max] = max(psd_mild_mean);
    F_maxPSD = round(F_AOI(idx_max));
    plot([F_maxPSD F_maxPSD], [0 maxPSD + maxPSD * 0.2], 'k--')

    xlim([min(F_AOI) max(F_AOI)])

    % legend
    legend([h1, h2, h3], {'normal', 'mild', 'moderate'})

    % title
    title(['PSD in ' brainarea])

    % save figure
    savename = fullfile(savefolder, ['psd_' brainarea]);
    saveas(gcf, savename, 'png')

    clear psd_allsegs_normal psd_allsegs_mild psd_allsegs_moderate
    clear psd_normal_FAOI psd_mild_FAOI psd_moderate_FAOI
    clear psd_normal_high psd_normal_low psd_normal_mean psd_mild_high psd_mild_low psd_mild_mean psd_moderate_high psd_moderate_low psd_moderate_mean
    clear h1 h2 h3 maxPSD F_maxPSD idx_max
    clear savename

end

function plotPSD_comp_multichns(psd_normal, psd_mild, psd_moderate, F_all, plotF_AOI, savefolder, brainarea)
    %%  plot the psd comparison of normal and mild
    %
    %   Inputs:
    %       psd_normal, psd_mild, psd_moderate: psd of all segments in normal, mild and moderate, nfs * nchns * nsegs
    %
    %       F_all: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
    %
    %       plotF_AOI: the plot frequence range (in hertz)  (1 * 2)
    %
    %       savefile_prefix: the savefile prefix including folder

    % find the idx for F_AOI
    idx_AOI = find(F_all >= plotF_AOI(1) & F_all <= plotF_AOI(2));

    % colors setup
    color_normal_range = [224, 255, 255] / 255;
    color_normal_mean = [0, 0, 255] / 255;
    color_mild_range = [255, 228, 225] / 255;
    color_mild_mean = [255, 00, 0] / 255;
    color_moderate_range = [238, 238, 238] / 255;
    color_moderate_mean = [0, 0, 0] / 255;

    % plot setup
    linewidth = 1.5;

    % F_AOI
    F_AOI = F_all(idx_AOI);
    % reshape F_AOI
    [m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m * n); clear m n

    nchns = size(psd_normal, 2);

    for chni = 1:nchns

        psd_allsegs_normal = squeeze(psd_normal(:, chni, :));
        psd_allsegs_mild = squeeze(psd_mild(:, chni, :));
        psd_allsegs_moderate = squeeze(psd_moderate(:, chni, :));

        % extract the  psd of AOI
        psd_normal_FAOI = psd_allsegs_normal(idx_AOI, :);
        psd_mild_FAOI = psd_allsegs_mild(idx_AOI, :);
        psd_moderate_FAOI = psd_allsegs_moderate(idx_AOI, :);

        psd_normal_high = max(psd_normal_FAOI, [], 2);
        psd_normal_low = min(psd_normal_FAOI, [], 2);
        psd_normal_mean = mean(psd_normal_FAOI, 2);

        psd_mild_high = max(psd_mild_FAOI, [], 2);
        psd_mild_low = min(psd_mild_FAOI, [], 2);
        psd_mild_mean = mean(psd_mild_FAOI, 2);

        psd_moderate_high = max(psd_moderate_FAOI, [], 2);
        psd_moderate_low = min(psd_moderate_FAOI, [], 2);
        psd_moderate_mean = mean(psd_moderate_FAOI, 2);

        % reshape, psd_*_high/low/mean into 1 * nfs

        [m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m * n); clear m n
        [m, n] = size(psd_normal_low); psd_normal_low = reshape(psd_normal_low, 1, m * n); clear m n
        [m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m * n); clear m n

        [m, n] = size(psd_mild_high); psd_mild_high = reshape(psd_mild_high, 1, m * n); clear m n
        [m, n] = size(psd_mild_low); psd_mild_low = reshape(psd_mild_low, 1, m * n); clear m n
        [m, n] = size(psd_mild_mean); psd_mild_mean = reshape(psd_mild_mean, 1, m * n); clear m n

        [m, n] = size(psd_moderate_high); psd_moderate_high = reshape(psd_moderate_high, 1, m * n); clear m n
        [m, n] = size(psd_moderate_low); psd_moderate_low = reshape(psd_moderate_low, 1, m * n); clear m n
        [m, n] = size(psd_moderate_mean); psd_moderate_mean = reshape(psd_moderate_mean, 1, m * n); clear m n


        % plot range
        figure
        set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
        fill([F_AOI flip(F_AOI)], [psd_normal_high flip(psd_normal_low)], color_normal_range, 'LineStyle', 'none')
        hold all
        fill([F_AOI flip(F_AOI)], [psd_mild_high flip(psd_mild_low)], color_mild_range, 'LineStyle', 'none')
        fill([F_AOI flip(F_AOI)], [psd_moderate_high flip(psd_moderate_low)], color_moderate_range, 'LineStyle', 'none')

        % plot mean
        h1 = plot(F_AOI, psd_normal_mean, 'Color', color_normal_mean, 'LineWidth', linewidth);
        h2 = plot(F_AOI, psd_mild_mean, 'Color', color_mild_mean, 'LineWidth', linewidth);
        h3 = plot(F_AOI, psd_moderate_mean, 'Color', color_moderate_mean, 'LineWidth', linewidth);

        % find the frequency with maximum density
        [maxPSD, idx_max] = max(psd_mild_mean);
        F_maxPSD = round(F_AOI(idx_max));
        plot([F_maxPSD F_maxPSD], [0 maxPSD + maxPSD * 0.2], 'k--')

        xlim([min(F_AOI) max(F_AOI)])

        % legend
        legend([h1, h2, h3], {'normal', 'mild', 'moderate'})

        % title
        title(['PSD  in ' upper(brainarea) num2str(chni-1) '-' num2str(chni)])

        % save figure
        savename = fullfile(savefolder, ['psd_' brainarea '_ch' num2str(chni)]);
        saveas(gcf, savename, 'png')

        clear psd_allsegs_normal psd_allsegs_mild psd_allsegs_moderate 
        clear psd_normal_FAOI psd_mild_FAOI psd_moderate_FAOI
        clear psd_normal_high psd_normal_low psd_normal_mean psd_mild_high psd_mild_low psd_mild_mean psd_moderate_high psd_moderate_low psd_moderate_mean
        clear h1 h2 h3  maxPSD F_maxPSD idx_max
        clear savename
    end

end
