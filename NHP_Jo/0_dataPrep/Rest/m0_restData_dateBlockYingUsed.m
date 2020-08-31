function m0_restData_dateBlockYingUsed()
    %   Manually marked the good and bad segments
    %
    % 1. add variable segsRemain, for marking each segment with 1 (good) or 0 (not good)
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
    [codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);
    [datafolder, ~,~,~] = exp_subfolders();


    %% global parameter
    animal = 'Jo';

    %% save setup
    savefilename = 'dateBlocksYingUsed_rest.xlsx';
    savefolder = codecorresfolder;
    savefile = fullfile(savefolder, savefilename);


    if exist(savefile, 'file') % savefile exist
        disp([savefile 'exist!'])
        reply = input('Do you want to rerun (load large data_segments file)? y/n [n]:', 's');
        if isempty(reply)
            reply = 'n';
        end
        reply = lower(reply);
        if reply == 'n'
            return;         
        end
    end


    %% Start Here
    file_datasegs_Ying = 'data_segment_Jo_Rest_NonMovwithMAbyCleandataNo60HzFilt (2).mat';
    folder_local = fullfile(datafolder, animal); % folder_local = '/home/lingling/Desktop/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/data/Jo'

    if exist(fullfile(folder_local,  file_datasegs_Ying), 'file')
        disp('Loading data_segments file locally.....')
        load(fullfile(folder_local, file_datasegs_Ying),'data_segments')
    else
        disp('Loading data_segments from server.......')
        load(fullfile('/home/lingling/root2/Ying Yu/SKB_EventRelatedBeta/datasets_Allanimals/Resting/', file_datasegs_Ying), 'data_segments');
    end


    nfiles = length(data_segments);
    for i = 1:nfiles

        dateBlock = getfield(data_segments, {i}, 'date');
        dbs_info = getfield(data_segments, {i}, 'dbs_info');

        if (~strcmp(dbs_info, 'n/a')) % only extract the dateblocks without any dbs stimuli
            continue;
        end

        % extract condition
        tmp = strsplit(dateBlock, '_');
        dateofexp = datenum(tmp{1}, 'yyyymmdd');
        condition = parsePDCondition(dateofexp, animal);
        clear tmp dateofexp

        if i == 1
            dateBlocks = dateBlock;
            predate = dateBlock;
            cond = {condition};
        else

            if (~strcmp(predate, dateBlock))
                dateBlocks = cat(1, dateBlocks, dateBlock);
                predate = dateBlock;
                cond = [cond; condition];
            end

        end

        clear dateBlock
    end


    tbl_dateblocks = table(dateBlocks, cond, 'VariableNames', {'dateBlock_rest', 'cond'});

    writetable(tbl_dateblocks, savefile);
    disp(['Extracted dateblockstr in ' savefile])
    clear tbl_dateblocks dateBlock_rest cond
end
