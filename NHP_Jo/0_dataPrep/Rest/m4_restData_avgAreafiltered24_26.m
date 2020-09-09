function m4_restData_avgAreafiltered24_26()
    %% narrow filtered rest data recorded with grayMatter in frequency [26 28]Hz
    %
    %   1. averaged the lfp in one area (except the DBS)
    %
    %   2. narrowed filltered.
    %
    %   3. seg into intervals with same time length
    %
    %
    %
    %   Input:
    %
    %       \m2_restData_selectSeg_Power
    %
    %       'data_segments', 'fs', 'T_chnarea'
    %
    %
    %   Output variables:
    %
    %
    %        lfpsegs: lfp segments in all areas
    %
    %        fs: resampled samping rate (default 500Hz)
    %
    %        chnAreas: chnAreas cell for used in python

    %% folders generate
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');

    % find the codefolder
    idx = strfind(codefilepath, 'code');
    codefolder = codefilepath(1:idx + length('code') - 1);
    clear idx

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));

    % the corresponding pipeline and the parent folder for this code
    [codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

    %% global variables
    % animal
    tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
    animal = tmp(length('/NHP_') + 1:end - 1);

    segt = 2;

    %%  input setup
    % band pass frequency
    frebp = [24 26];
    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_Power');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = ['filtered' num2str(frebp(1)) '_' num2str(frebp(2)) '_seg_avgArea'];

    %% starting: narrow filter the lfp data of all the files
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);
    f = waitbar(0, ['Narrow Filtering....']);

    for filei = 1:nfiles
        % wait bar
        waitbar(filei / nfiles, f, ['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);

        filename = files(filei).name;

        % extract filtered lfp with same length from 1 file
        [lfpsegs, T_chnsarea, fs] = avgAreafilter_seg_1file(fullfile(inputfolder, filename), frebp, segt);

        % extract chnAreas cell for used in python
        chnAreas = T_chnsarea.brainarea;

        idx_STN = find(strcmp(T_chnsarea.brainarea, 'STN'));
        for i = 1:length(idx_STN)
            chnAreas{idx_STN(i)} = ['stn' num2str(i - 1) '-' num2str(i)];
        end

        idx_GP = find(strcmp(T_chnsarea.brainarea, 'GP'));
        for i = 1:length(idx_GP)
            chnAreas{idx_GP(i)} = ['gp' num2str(i - 1) '-' num2str(i)];
        end

        % save
        idx = strfind(filename, [animal '_']);
        tmpn = length([animal '_']);
        savefilename = [filename(idx:idx + tmpn - 1) savefilename_addstr ...
                        upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];

        save(fullfile(savefolder, savefilename), 'lfpsegs', 'chnAreas', 'fs');

        clear filename lfpsegs T_chnsarea fs
        clear idx tmpn savefilename

    end

    close(f);
    disp(['narrow filtered lfpdata with same length are saved to ' savefolder])
end

function [filted_lfp_segs, T_chnsarea_new, fs] = avgAreafilter_seg_1file(file, frebp, segt)
    %
    % Args:
    %   file:  one file (full path)
    %   frebp: bandpass filter range (e.g.[5 50])
    %   segt: the same time length, unit second (e.g 2)
    %
    %  Returns:
    %   filted_lfp_segs: bandpass filted lfp with same length (ntemp * nchns * nsegs)
    %   T_chnsarea_new: new T_chnsarea_new (height = nchns)
    %   fs: sample rate

    load(file, 'data_segments', 'fs', 'T_chnsarea')

    uniqBrainAreas = unique(T_chnsarea.brainarea);

    ntemp_same = ceil(fs * segt);
    filted_lfp_segs = []; % filted_lfp_segs: ntemp * nchns * nsegs
    T_chnsarea_new = T_chnsarea([], :);

    for segi = 1:length(data_segments)

        %%% band pass filter, for not DBS, average first%%%
        filted_lfp_1seg = [];

        for areai = 1:length(uniqBrainAreas)
            brainarea = uniqBrainAreas{areai};

            % extract the idx of brainarea
            mask_area = strcmp(T_chnsarea.brainarea, brainarea);

            % extract the lfp data of segi in brainarea, lfp_oneseg: ntemp * nchns_area
            lfp_oneseg = data_segments(segi).lfp(:, mask_area);

            if ~strcmp(brainarea, 'STN') &&~strcmp(brainarea, 'GP')% if not DBS case, use averaged lfp

                lfp_oneseg = mean(lfp_oneseg, 2);

            end

            %%% band pass filter, filted_lfp_area: ntemp * 1/nchns_STN %%%
            filted_lfp_area = bp_filter(lfp_oneseg, frebp, fs);

            % concatenate along chns
            filted_lfp_1seg = cat(2, filted_lfp_1seg, filted_lfp_area);

            %%% new T_chnsarea%%%%
            if (segi == 1)% only generate new T_chnsarea once

                T_area = T_chnsarea(mask_area, :);

                if ~strcmp(brainarea, 'STN') &&~strcmp(brainarea, 'GP')% only has one row if not DBS
                    T_area = T_area(1, :);
                    T_area.recordingchn = NaN;
                end

                T_chnsarea_new = [T_chnsarea_new; T_area];
                clear T_area
            end

            clear brainarea mask_area T_area lfp_oneseg filted_lfp_area

        end

        %%% seg into intervals with same time length   %%%%
        for segi_same = 1:floor(size(filted_lfp_1seg, 1) / ntemp_same)
            idx_str = (segi_same - 1) * ntemp_same + 1;
            idx_end = segi_same * ntemp_same;

            filted_lfp_segs = cat(3, filted_lfp_segs, filted_lfp_1seg(idx_str:idx_end, :));

            clear idx_str idx_end
        end

        clear filted_lfp_1seg segi_same areai

    end

    T_chnsarea_new.chni = [1:height(T_chnsarea_new)]';
end

function filteredX = bp_filter(X, frebp, fs)
    %% band pass for each channel X: ntemporal * nchns

    [ntemp, nchns] = size(X);

    if ntemp == 1
        X = X';
        [ntemp, nchns] = size(X);
    end

    filteredX = zeros(ntemp, nchns);

    for chi = 1:nchns
        filteredX(:, chi) = filter_bpbutter(X(:, chi), frebp, fs);
    end

end
