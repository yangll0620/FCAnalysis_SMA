function m2_SKTData_usefulChannel()


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup

% input folder: extracted raw STK data
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_preprocessing');


% depth and reference channel for M1
depth_M1Layer5 = [1.25 1.75] * 8;
recordchnref_M1 = 58;

% depth and reference channel for PMC
depth_PMCLayer5 = [1.25 1.75] * 8;
recordchnref_PMC = 50;


%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'usefulChannels';




%% Code Start Here
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
close all;
f = waitbar(0, ['Processing lfp data...']);
for i = 1 : nfiles
    % wait bar
    waitbar(i/nfiles,f,['Processing  file ' num2str(i) '/' num2str(nfiles)]);
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(i).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    disp(filename)
    
    
    % extract the date of the exp
    idx = strfind(filename, '_bktdt');
    dateofexp = datetime(filename(idx-6: idx-1), 'InputFormat', 'MMddyy');
    
    
    %%% ----- extract useful lfpm1 and recordchnUsed_M1  ------%%%
    brainarea = 'M1';
    depth_useful  = depth_M1Layer5;
    recordchnref = recordchnref_M1;
    
    
    % extract recording number for brainarea in Useful Layer on dateofexp
    recordchnUsed_dateofexp = useful_recordchn_extract(brainarea, depth_useful, dateofexp);
    if isempty(recordchnUsed_dateofexp)
        disp(["ignore this as useful channel of " + brainarea + " is empty"])
        continue;
    end
    
    % extract the reference lfp in brainarea
    chansRef = ismember(T_chnsarea.recordingchn, recordchnref);
    lfp_ref = lfpdata(chansRef, :,:);
    
    % extract the bipolar lfp for each channel in Useful Layer
    chansUsed = ismember(T_chnsarea.recordingchn, recordchnUsed_dateofexp);
    lfp_area = lfpdata(chansUsed, :,:)  - repmat(lfp_ref, length(find(chansUsed)), 1, 1);
    
    % reassign
    lfpM1 = mean(lfp_area, 1);
    recordchnUsed_M1 = recordchnUsed_dateofexp;
    clear chansRef lfp_area lfp_ref chansUsed recordchnUsed_dateofexp brainarea
    
    
    
    
    %%% ----- extract useful lfpPMC and recordchnUsed_PMC  ------%%%
    brainarea = 'PMC';
    depth_useful  = depth_PMCLayer5;
    recordchnref = recordchnref_PMC;
    
    
    % extract recording number for brainarea in Useful Layer on dateofexp
    recordchnUsed_dateofexp = useful_recordchn_extract(brainarea, depth_useful, dateofexp);
    if isempty(recordchnUsed_dateofexp)
        disp(["ignore this as useful channel of " + brainarea + " is empty"])
        continue;
    end
    
    % extract the reference lfp in brainarea
    chansRef = ismember(T_chnsarea.recordingchn, recordchnref);
    lfp_ref = lfpdata(chansRef, :,:);
    
    % extract the bipolar lfp for each channel in Useful Layer
    chansUsed = ismember(T_chnsarea.recordingchn, recordchnUsed_dateofexp);
    lfp_area = lfpdata(chansUsed, :,:)  - repmat(lfp_ref, length(find(chansUsed)), 1, 1);
    
    % reassign
    lfpPMC = mean(lfp_area, 1);
    recordchnUsed_PMC = recordchnUsed_dateofexp;
    clear chansRef lfp_area lfp_ref chansUsed recordchnUsed_dateofexp brainarea
    
    
    
    
    
    %%% ----- deal with lfpdata and T_chnsarea  ------%%%
    rowidx_M1PMC = (T_chnsarea.brainarea == "M1" | T_chnsarea.brainarea == "PMC");
    lfpdata(rowidx_M1PMC, :, :) = [];
    lfpdata = cat(1, lfpM1, lfpPMC, lfpdata);
    
    %%%%%%%  T_chnsarea 
    % delete the rows for M1 and PMC
    T_chnsarea(rowidx_M1PMC,:) = [];
    
    % generate two new rows for M1 and PMC
    newrows_M1PMC = T_chnsarea(1:2,:);
    % recordingchn to be nan
    newrows_M1PMC.recordingchn(1:2) = nan;
    % record chn Used stored in notes
    newrows_M1PMC.notes{1} = recordchnUsed_M1;
    newrows_M1PMC.notes{2} = recordchnUsed_PMC;
    % record chn Used stored in notes
    newrows_M1PMC.brainarea{1} = 'M1';
    newrows_M1PMC.brainarea{2} = 'PMC';
    
    
    % add the two new rows
    T_chnsarea = [newrows_M1PMC;T_chnsarea];
    T_chnsarea.chni = [1:height(T_chnsarea)]';
    
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    clear newrows_M1PMC
    clear lfpM1 lfpPMC recordchnUsed_PMC recordchnUsed_M1 rowidx_M1PMC
    clear lfpdata fs T_chnsarea T_idxevent
    
end

close(f)
disp("Processing all files Done!")





