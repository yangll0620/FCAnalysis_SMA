function m3_restData_PSDDBS_extract()
    %   PSD estimates for mild and normal individually
    %
    %   psd for each brain area, as well as each DBS contact

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

    % pwelch psd estimate variable
    twin_pwelch = 2;

    % variables for plotting
    plotF_AOI = [5 50];

    % input folder
    inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_M1Power');

    %% save setup
    savefolder = codecorresfolder;
    savefilename_prefix = 'psd_';


    

end
