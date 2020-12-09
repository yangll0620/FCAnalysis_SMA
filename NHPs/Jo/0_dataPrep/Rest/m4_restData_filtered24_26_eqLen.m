function m4_restData_filtered24_26_eqLen()
    %% narrow filtered rest data recorded with grayMatter in frequency [26 28]Hz
    %
    %   1. narrowed filltered.
    %
    %   2. seg into intervals with same time length
    %
    %   3. ignore files with no segment
    %
    %
    %   Input:
    %
    %       \m3_restData_rmChns_avgArea
    %
    %       'data_segments', 'fs', 'T_chnarea'
    %
    %
    %   Output variables:
    %
    %
    %        lfpsegs: lfp segments in all areas
    %
    %        fs: samping rate (default 500Hz)
    %
    %        chnAreas: chnAreas cell for used in python
    %

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
    [fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
    animal = codecorresfolder(fi + length('NHPs/'):j);

    segt = 2;

    %%  input setup
    % band pass frequency
    frebp = [24 26];
    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm3_restData_rmChns_avgArea');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = ['filtered' num2str(frebp(1)) '_' num2str(frebp(2)) '_eqLen'];

    %% starting: narrow filter the lfp data of all the files
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);
    f = waitbar(0, ['Narrow Filtering....']);

    for filei = 1:nfiles
        % wait bar
        waitbar(filei / nfiles, f, ['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);

        filename = files(filei).name;

        % extract filtered lfp with same length from 1 file
        [lfpsegs, T_chnsarea, fs] = filter_seg_1file(fullfile(inputfolder, filename), frebp, segt);

        if isempty(lfpsegs)
            continue;
        end

        % extract chnAreas cell for used in python
        chnAreas = T_chnsarea.brainarea;


        % save
        [i, j]= regexp(filename, [animal '_[a-zA-Z]*_' ]);
        savefilename = [animal '_' savefilename_addstr filename(j:end)];

        save(fullfile(savefolder, savefilename), 'lfpsegs', 'chnAreas', 'fs');

        clear filename lfpsegs T_chnsarea fs
        clear idx tmpn savefilename

    end

    close(f);
    disp(['narrow filtered lfpdata with same length are saved to ' savefolder])
end

function [filted_eqLen_lfp, T_chnsarea, fs] = filter_seg_1file(file, frebp, segt)
    %
    % Args:
    %   file:  one file (full path)
    %   frebp: bandpass filter range (e.g.[5 50])
    %   segt: the same time length, unit second (e.g 2)
    %
    %  Returns:
    %   filted_eqLen_lfp: bandpass filted lfp with same length (nchns * ntemp  * nsegs)
    %   T_chnsarea
    %   fs: sample rate

    disp(file)
    load(file, 'data_segments', 'fs', 'T_chnsarea')

    % extract uniqBrainAreas
    mask_emptyarea = cellfun(@(x) isempty(x), T_chnsarea.brainarea);
    uniqBrainAreas = unique(T_chnsarea.brainarea(~mask_emptyarea));

    ntemp_same = ceil(fs * segt);
    filted_eqLen_lfp = []; % filted_eqLen_lfp: nchns * ntemp * nsegs


    for segi = 1:length(data_segments)

        %%% band pass filter, for not DBS, average first%%%
        filted_lfp_1seg = [];

        for areai = 1:length(uniqBrainAreas)
            brainarea = uniqBrainAreas{areai};

            % extract the idx of brainarea
            mask_area = strcmp(T_chnsarea.brainarea, brainarea);

            % extract the lfp data of segi in brainarea, lfp_oneseg: nchns_area * ntemp
            lfp_oneseg = data_segments(segi).lfp(mask_area, :);


            %%% band pass filter, filted_lfp_area: ntemp * 1/nchns_STN %%%
            filted_lfp_area = bp_filter(lfp_oneseg, frebp, fs);

            % concatenate along chns
            filted_lfp_1seg = cat(1, filted_lfp_1seg, filted_lfp_area);

            clear brainarea mask_area T_area lfp_oneseg filted_lfp_area

        end

        %%% seg into intervals with same time length   %%%%
        for segi_same = 1:floor(size(filted_lfp_1seg, 2) / ntemp_same)
            idx_str = (segi_same - 1) * ntemp_same + 1;
            idx_end = segi_same * ntemp_same;

            filted_eqLen_lfp = cat(3, filted_eqLen_lfp, filted_lfp_1seg(:, idx_str:idx_end));

            clear idx_str idx_end
        end

        clear filted_lfp_1seg segi_same areai

    end

end

function filteredX = bp_filter(X, frebp, fs)
    %% band pass for each channel X: nchns *  ntemp

    [nchns, ntemp] = size(X);

    if ntemp == 1
        X = X';
        [nchns, ntemp] = size(X);
    end

    filteredX = zeros(nchns, ntemp);

    for chi = 1:nchns
        filteredX(chi, :) = filter_bpbutter(X(chi, :), frebp, fs);
    end

end
