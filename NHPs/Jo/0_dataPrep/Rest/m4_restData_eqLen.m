function m4_restData_eqLen()
    %%  Seg into equal length 
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

    % input folder: extracted raw rest data with grayMatter
    inputfolder = fullfile(codecorresParentfolder, 'm3_restData_rmChns_avgArea');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_addstr = 'eqLen';

    %% starting: narrow filter the lfp data of all the files
    files = dir(fullfile(inputfolder, '*.mat'));
    nfiles = length(files);
    
    for filei = 1:nfiles
        
        filename = files(filei).name;
        
        load(fullfile(inputfolder, filename), 'data_segments', 'fs', 'T_chnsarea')
        
        %%% seg into equal length %%% 
        lfpsegs = [];
        ntemp_same = ceil(fs * segt);
        for segi = 1:length(data_segments)
            lfp_1seg = data_segments(segi).lfp;
            %%% seg into intervals with same time length   %%%%
            for segi_same = 1:floor(size(lfp_1seg, 2) / ntemp_same)
                idx_str = (segi_same - 1) * ntemp_same + 1;
                idx_end = segi_same * ntemp_same;
                lfpsegs = cat(3, lfpsegs, lfp_1seg(:, idx_str:idx_end));
                clear idx_str idx_end
            end
        end
        clear ntemp_same
        
        
        % extract chnAreas cell for used in python
        chnAreas = T_chnsarea.brainarea;
        
        
        % save
        [~, j]= regexp(filename, [animal '_[a-zA-Z]*_' ]);
        savefilename = [animal '_' savefilename_addstr filename(j:end)];

        save(fullfile(savefolder, savefilename), 'lfpsegs', 'chnAreas', 'fs');
        
        clear lfpsegs chnAreas fs T_chnsarea
        clear filename
    end
