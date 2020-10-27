function m4_restData_filtered13_15()
    %% narrow filtered rest data recorded with grayMatter in frequency [13 15]Hz
    %
    %
    %   narrowed filltered.
    %
    %
    %   Input:
    %
    %       m3_restData_eqLen_avgArea_combVLoVPLo
    %
    %       'lfpdata', 'fs', 'T_chnarea'
    %
    %
    %   Output variables:
    %
    %
    %        lfpdata: filtered lfp data in all areas
    %
    %        fs: sample rate
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
    

    %%  input setup
    % band pass frequency
    frebp = [13 15];
    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm3_restData_eqLen_avgArea_combVLoVPLo');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = ['filtered' num2str(frebp(1)) '_' num2str(frebp(2))];

    %% starting: narrow filter the lfp data of all the files
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);
    f = waitbar(0, ['Narrow Filtering....']);

    for filei = 7:nfiles
        % wait bar
        waitbar(filei / nfiles, f, ['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);

        filename = files(filei).name;
        
        % lfpdata: nareas * ntemp * nsegs
        load(fullfile(files(filei).folder, filename), 'lfpdata', 'fs','T_chnsarea');
        
        
        %%% narrow filtered %%%
        filtered_lfpdata = bp_filter(lfpdata, frebp, fs);
        

        %%% extract chnAreas cell for used in python %%%
        chnAreas = T_chnsarea.brainarea;

        % change STN into stn0-1, stn1-2 et.al
        idx_STN = find(strcmp(T_chnsarea.brainarea, 'STN'));
        for i = 1:length(idx_STN)
            chnAreas{idx_STN(i)} = ['stn' num2str(i - 1) '-' num2str(i)];
        end

        % change GP into gp0-1, gp1-2 et.al
        idx_GP = find(strcmp(T_chnsarea.brainarea, 'GP'));
        for i = 1:length(idx_GP)
            chnAreas{idx_GP(i)} = ['gp' num2str(i - 1) '-' num2str(i)];
        end
        
        
        %%% save  %%%
        lfpdata = filtered_lfpdata;
        
        idx = [strfind(filename, '_normal_') strfind(filename, '_mild_') strfind(filename, '_moderate_')];
        savefilename = [animal '_' savefilename_addstr filename(idx:end)];

        save(fullfile(savefolder, savefilename), 'lfpdata', 'chnAreas', 'fs');

        clear filename lfpdata T_chnsarea fs chnAreas
        clear filtered_lfpdata  idx_STN idx_GP
        clear idx tmpn savefilename

    end

    close(f);
    disp(['narrow filtered lfpdata are saved to ' savefolder])
end


function filteredX = bp_filter(X, frebp, fs)
    %% band pass for each channel X: nareas * ntemp * nsegs

    [nareas, ntemp, nsegs] = size(X);

    filteredX = zeros(nareas, ntemp, nsegs);
    for ai = 1:nareas
        for segi = 1: nsegs
            filteredX(ai, :, segi) = filter_bpbutter(X(ai, :, segi), frebp, fs);
        end
    end

end
