function m0_restData_dateBlocksYingUsed()

    %% folder generate

    % the full path and the name of code file without suffix
    codefilepath = mfilename('fullpath');

    % find the codefolder
    idx = strfind(codefilepath, 'code');
    codefolder = codefilepath(1:idx + length('code') - 1);
    clear idx

    % add util path
    addpath(genpath(fullfile(codefolder, 'util')));

    %% global parameter
    animal = 'Pinky';

    % the corresponding pipeline folder for this code
    [codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

    % datafolder
    [datafolder, ~, ~, ~] = exp_subfolders();

    %% input setup
    % dates used by Ying
    filesYingUsed = fullfile(datafolder, 'Pinky', 'DateYingUsed_Pinky_ForYLL.mat');

    %% save setup
    savefolder = codecorresfolder;
    savefilename = [animal 'dateBlocksYingUsed_rest.xlsx'];
    savefile = fullfile(savefolder, savefilename);

    %% Start Here

    load(filesYingUsed, 'dateYingAnalyzed')

    for i = 1:size(dateYingAnalyzed, 1)
        dateBlock = string(dateYingAnalyzed(i, :));

        % extract condition
        tmp = strsplit(dateBlock, '_');
        dateofexp = datenum(tmp{1}, 'yyyymmdd');
        condition = parsePDCondition(dateofexp, animal);
        clear tmp dateofexp

        if i == 1
            dateBlocks = dateBlock;
            predate = dateBlock;
            conditions = {condition};
        else

            if (predate ~= dateBlock)
                dateBlocks = cat(1, dateBlocks, dateBlock);
                predate = dateBlock;
                conditions = [conditions; condition];
            end

        end

        clear dateBlock
    end

    tbl_dateblocks = table(dateBlocks, conditions, 'VariableNames', {'dateBlock_rest', 'condition'});

    writetable(tbl_dateblocks, savefile);
    disp(['Extracted dateblockstr in ' savefile])
    clear tbl_dateblocks dateBlock_rest conditions
end