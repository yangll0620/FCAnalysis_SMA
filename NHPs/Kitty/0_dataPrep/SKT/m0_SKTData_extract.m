function m0_SKTData_extract()    
    %% extract the COT normal data and the SKT moderate data for Kitty 

    
    
    
    %% extract the corresponding pipeline folder for this code
    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');
    % code folder
    codefolder = codefilepath(1:strfind(codefilepath, 'code') + length('code') - 1);

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));


    % datafolder, pipelinefolder
    [codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);
    
    
    
    %% combine normalCOT_folder and moderateSKT_folder into m0_SKTData_extract
    normalCOT_folder = fullfile(codecorresParentfolder, 'm0_normalCOTData_extract');
    moderateSKT_folder = fullfile(codecorresParentfolder, 'm0_moderateSKTData_extract');
    
    files = [dir(fullfile(normalCOT_folder, '*.mat')); dir(fullfile(moderateSKT_folder, '*.mat'))];
    for filei = 1 : length(files)
        src = fullfile(files(filei).folder, files(filei).name);
        des = fullfile(codecorresfolder, files(filei).name);
        
        copyfile(src, des);
        
        clear src des
    end
    
    
end