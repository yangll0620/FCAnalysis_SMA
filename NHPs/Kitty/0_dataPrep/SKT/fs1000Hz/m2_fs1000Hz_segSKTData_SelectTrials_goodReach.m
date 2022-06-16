function m2_fs1000Hz_segSKTData_SelectTrials_goodReach()
%  select trials for fs1000 using fs500 data
% 


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);



%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder_fs1000 = fullfile(codecorresParentfolder, 'm1_fs1000Hz_segSKTData_avgArea');


inputfolder_fs500 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV';
inputfolder_fs500_avgArea = fullfile(inputfolder_fs500, 'm1_segSKTData_avgArea');
inputfolder_fs500_selectedTrials = fullfile(inputfolder_fs500, 'm2_segSKTData_SelectTrials_goodReach');

%% save setup
savefolder = codecorresfolder;

%% code start here
files = dir(fullfile(inputfolder_fs1000, '*.mat'));
for fi = 1 : length(files)
    filename = files(fi).name;
    data_fs500 = load(fullfile(inputfolder_fs500_avgArea, filename));
    load(fullfile(inputfolder_fs1000, filename), 'T_chnsarea', 'T_idxevent_ma', 'Wrist_smooth_trial', 'fs_ma', 'mask_goodreach', 'smoothWspeed_trial',...  
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
    
    selectedFile = dir(fullfile(inputfolder_fs500_selectedTrials, ['*' filename(end-19:end)]));
    if(length(selectedFile)~=1)
        disp([filename ': selected File not 1.']);
        clear selectedFile
        continue
    end
    
    load(fullfile(inputfolder_fs500_selectedTrials, selectedFile(1).name), 'idxGroups', 'tbl_selectedTrialsMarks', 'selectedTrials');
    save(fullfile(savefolder, selectedFile(1).name), 'idxGroups', 'tbl_selectedTrialsMarks', 'selectedTrials', ...
                                                                       'T_chnsarea', 'T_idxevent_ma', 'Wrist_smooth_trial', 'fs_ma', 'mask_goodreach', 'smoothWspeed_trial',...  
                                                                       'T_idxevent_lfp', 'Wpos_smooth_trial', 'fs_lfp', 'lfpdata','mask_goodreturn');
    
    clear clear filename
    clear('idxGroups', 'tbl_selectedTrialsMarks', 'selectedTrials', ...
                                                                       'T_chnsarea', 'T_idxevent_ma', 'Wrist_smooth_trial', 'fs_ma', 'mask_goodreach', 'smoothWspeed_trial',...  
                                                                       'T_idxevent_lfp', 'Wpos_smooth_trial', 'fs_lfp', 'lfpdata','mask_goodreturn');
end


