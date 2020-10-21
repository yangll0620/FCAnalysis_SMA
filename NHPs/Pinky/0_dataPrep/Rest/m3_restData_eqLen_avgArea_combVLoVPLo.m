function m3_restData_eqLen_avgArea_combVLoVPLo()
    %% segment into intervals with same length (segt = 2 ) and average lfp across each area (combine VLo and VPLo)
    %
    %   remove unwanted chns lCd and rMC, and combine VLo and VPLo
    %
    % Input:
    %   m2_restData_selectSeg_Power
    %
    % Outputs:
    %   lfpdata: averaged lfp with same length ((nareas + nDBS) * ntemp * nsegs)
    %   T_chnsarea: new T_chnsarea_new (height = (nareas + nDBS))
    %   fs: sample rate 


    %% folders generate
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');

    % find the codefolder
    idx = strfind(codefilepath, 'code');
    codefolder = codefilepath(1:idx + length('code') - 1);
    clear idx

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));
    
    addpath(genpath(fullfile('.', 'subFuncs')));

    % the corresponding pipeline and the parent folder for this code
    [codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

    %% global variables
    % animal
    tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
    animal = tmp(length('/NHP_') + 1:end - 1);

    segt = 2;

    %%  input setup

    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_Power');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = 'eqLen_avgArea';

    %% starting: narrow filter the lfp data of all the files
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);
    
    
    for filei = 1:nfiles
        
        filename = files(filei).name;
        file = fullfile(files(filei).folder, files(filei).name);
           
        
        load(file, 'data_segments', 'fs', 'T_chnsarea')
        
        if isempty(data_segments)
            continue;
        end
        
        
        
        % replace *VLo and *VPLo with VLoVPLo
        mask_VLoVPLo = strcmp(T_chnsarea.brainarea, 'lVLo') | strcmp(T_chnsarea.brainarea, 'lVPLo');
        T_chnsarea{mask_VLoVPLo, 'brainarea'} = {'lVLo/VPLo'};
        mask_VLoVPLo = strcmp(T_chnsarea.brainarea, 'rVLo') | strcmp(T_chnsarea.brainarea, 'rVPLo');
        T_chnsarea{mask_VLoVPLo, 'brainarea'} = {'rVLo/VPLo'};
        
        
        disp(filename)
        [lfpdata, T_chnsarea_new] = eqLen_avgArea_1file(data_segments, T_chnsarea, segt, fs);
        

        
        
        % save
        T_chnsarea = T_chnsarea_new;
        
        idx = strfind(filename, [animal '_']);
        tmpn = length([animal '_']);
        savefilename = [filename(idx:idx + tmpn - 1) savefilename_addstr ...
            upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];
        
        save(fullfile(savefolder, savefilename), 'lfpdata',  'T_chnsarea', 'fs');
        
        
        clear filename file lfpdata T_chnsarea fs idx tempn savefilename
    end
end


