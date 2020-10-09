function m1_restData_cleaned_tmpmoderate_extract()
    % temporally use t_start and t_end through looking at videos by lingling to generate moderate segments.
    %
    %   t_seg stored in 'moderate_temp.csv'

    addpath('/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code/toolbox/NexMatlabFiles')

    clear

    % the t_start and t_end for good resting segments stored in .csv file, through looking at video by lingling
    T_tmpmoderate = readtable('moderate_temp.csv', 'Delimiter',',');

    % -- save setup -- %
    savefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/NHP_Bug/0_dataPrep/Rest/m1_restData_cleaned_extract';
    savefilename_prefix = 'Bug_cleanedRestData_moderate_';
    dateformat_save = 'yyyymmdd';

    for i = 1:height(T_tmpmoderate)

        % extract each record of date_num, tdtblock and ts_seg
        date_num = datenum(string(T_tmpmoderate.date(i)), 'yyyymmdd');
        tdtblock = T_tmpmoderate.tdtblock(i);
        ts_seg = str2num(T_tmpmoderate.ts_seg{i});


        %  extract data_segments and fs for each paire of date_num and tdtblock
        [data_segments, fs] = datasegments_1day(date_num, tdtblock, ts_seg);

        % --- save part  ----%
        savefile = fullfile(savefolder, [savefilename_prefix datestr(date_num, dateformat_save) '_tdt' num2str(tdtblock) '.mat']);
        save(savefile, 'data_segments', 'fs');

        clear date_num tdtblock ts_seg data_segments fs
    end

    copyfile('moderate_temp.csv', fullfile(savefolder, 'moderate_temp.csv'));

end

function [data_segments, fs] = datasegments_1day(date_num, tdtblock, ts_seg)
    %   return the data segments for a pair of date_num and tdtblock with time segment array ts_seg
    %

    oneday_folder = fullfile('/home/lingling/root2/Animals2/Bug/Recording/Processed/DataDatabase', ['Bug_' datestr(date_num, 'mmddyy')]);

    % --- extract lfp_array ----%
    lfp_folder = fullfile(oneday_folder, 'LFP', ['Block-' num2str(tdtblock)]);

    lfp_array = [];

    for chi = 1:96

        filename_pattern = fullfile(lfp_folder, ['*LFPch' num2str(chi) '.nex']);
        files = dir(filename_pattern);

        if length(files) ~= 1
            disp([filename_pattern ' has ' num2str(length(files)) ' file, skip !']);

            clear filename_pattern files
            break;
        end

        % read from nex file
        [nexData] = readNexFile(fullfile(files(1).folder, files(1).name));

        %
        lfp_struct = nexData.contvars{1, 1};

        if strcmp(lfp_struct.name, ['LFPch' num2str(chi)])% is LFPch1 struct

            % fs
            if ~exist('fs', 'var')
                fs = lfp_struct.ADFrequency;
            end

            if fs ~= lfp_struct.ADFrequency
                disp([files(1).name ', lfp_struct.ADFrequency != ' num2str(fs)]);
            end

            % combine the chn data
            lfp_array = cat(2, lfp_array, lfp_struct.data);
        else
            disp([files(1).name ' name = ' lfp_struct.name]);
        end

        clear nexData filename_pattern files lfp_struct

    end

    clear lfp_folder chi

    % --- extract lfp_dbs ----%
    dbs_folder = fullfile(oneday_folder, 'DBSLFP', ['Block-' num2str(tdtblock)]);
    dbsfile_pattern = fullfile(dbs_folder, '*DBSLFP.nex');

    files = dir(dbsfile_pattern);

    if length(files) ~= 1
        disp([dbsfile_pattern ' has ' num2str(length(files)) ' file, skip !']);

    end

    % load nex dbs file
    [nexData] = readNexFile(fullfile(files(1).folder, files(1).name));

    lfp_dbs = [];

    for i = 3:18
        dbs_struct = nexData.contvars{i, 1};

        if strcmp(dbs_struct.name, ['RAW_DBSch' num2str(i - 2)])% is DBSch1 struct

            % check fs
            if fs ~= dbs_struct.ADFrequency
                disp([files(1).name ', dbs_struct.ADFrequency != ' num2str(fs)]);
            end

            % concatenate dbs along channel
            lfp_dbs = cat(2, lfp_dbs, dbs_struct.data);

        else
            disp([files(1).name ' nexData.contvars{' num2str(i) ', 1}.name is not ' 'RAW_DBSch' num2str(i - 1)]);

            clear dbs_struct

            break;
        end

        clear dbs_struct
    end

    clear nexData i

    % --- extract data_segments based on ts_seg ----%
    for segi = 1:size(ts_seg, 1)
        idx_segStr = round(ts_seg(segi, 1) * fs);
        idx_segEnd = round(ts_seg(segi, 2) * fs);

        % seg from lfp_array
        if ~exist('data_segments')
            data_segments(1).lfp_array = lfp_array(idx_segStr:idx_segEnd, :);
        else
            data_segments(end + 1).lfp_array = lfp_array(idx_segStr:idx_segEnd, :);
        end

        % lfp_stn and lfp_gp from lfp_dbs
        data_segments(end).lfp_stn = lfp_dbs(idx_segStr:idx_segEnd, 1:8);
        data_segments(end).lfp_gp = lfp_dbs(idx_segStr:idx_segEnd, 9:16);

        clear idx_segStr idx_segEnd
    end

    clear segi lfp_array lfp_dbs
end
