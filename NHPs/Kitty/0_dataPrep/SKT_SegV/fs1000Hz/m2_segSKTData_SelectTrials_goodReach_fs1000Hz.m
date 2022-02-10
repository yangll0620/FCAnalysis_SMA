folder_fs1000 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\fs1000Hz\m1_segSKTData_avgArea_fs1000Hz';
folder_fs500 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\m1_segSKTData_avgArea';
folder_selectedTrials_fs1000 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\fs1000Hz\m2_segSKTData_SelectTrials_goodReach_fs1000Hz';
folder_selectedTrials_fs500 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\m2_segSKTData_SelectTrials_goodReach';

files = dir(fullfile(folder_fs1000, '*.mat'));
for fi = 1 : length(files)
    filename = files(fi).name;
    data_fs500 = load(fullfile(folder_fs500, filename));
    load(fullfile(folder_fs1000, filename), 'T_chnsarea', 'T_idxevent_ma', 'Wrist_smooth_trial', 'fs_ma', 'mask_goodreach', 'smoothWspeed_trial',...  
                                                        'T_idxevent_lfp', 'Wpos_smooth_trial', 'fs_lfp', 'lfpdata','mask_goodreturn');
    
    if ~isequal(data_fs500.mask_goodreach, mask_goodreach)
        disp([filename ': mask_goodreach is not equal'])
        clear filename data_fs1000 data_fs500
        continue;
    end
    
    if ~isequal(data_fs500.mask_goodreturn, mask_goodreturn)
        disp([filename ': mask_goodreach is not equal'])
        clear filename data_fs1000 data_fs500
        continue;
    end
    clear data_fs1000 data_fs500
    
    selectedFile = dir(fullfile(folder_selectedTrials_fs500, ['*' filename(end-19:end)]));
    if(length(selectedFile)~=1)
        disp([filename ': selected File not 1.']);
        clear selectedFile
        continue
    end
    
    load(fullfile(folder_selectedTrials_fs500, selectedFile(1).name), 'idxGroups', 'tbl_selectedTrialsMarks', 'selectedTrials');
    save(fullfile(folder_selectedTrials_fs1000, selectedFile(1).name), 'idxGroups', 'tbl_selectedTrialsMarks', 'selectedTrials', ...
                                                                       'T_chnsarea', 'T_idxevent_ma', 'Wrist_smooth_trial', 'fs_ma', 'mask_goodreach', 'smoothWspeed_trial',...  
                                                                       'T_idxevent_lfp', 'Wpos_smooth_trial', 'fs_lfp', 'lfpdata','mask_goodreturn');
    
    clear clear filename
    clear('idxGroups', 'tbl_selectedTrialsMarks', 'selectedTrials', ...
                                                                       'T_chnsarea', 'T_idxevent_ma', 'Wrist_smooth_trial', 'fs_ma', 'mask_goodreach', 'smoothWspeed_trial',...  
                                                                       'T_idxevent_lfp', 'Wpos_smooth_trial', 'fs_lfp', 'lfpdata','mask_goodreturn');
end


