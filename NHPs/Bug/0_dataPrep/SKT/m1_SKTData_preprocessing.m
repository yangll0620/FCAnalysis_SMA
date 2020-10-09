function m1_SKTData_preprocessing()
%% preprocessing skt data 
%
%   Processing steps as follows:
%       1. add a new depth variable using the information from Bug_channelDepth_normal(mild).csv to T_chnsarea
%   
%       2. only extract lfp data in M1, SMA, PMC, thalamus (VA, VLo, VPLo) areas
%
%       3. bipolar for STN and GP channels
%
%       4. ignore the dateofexp which can not be added a depth variable due
%           to there is no row or more than two rows for the exp date
%
%       5. Down sample trials into fs_new = 500
%
%   Output variable:
%       fs: downsample rate fs = 500
%   
%       T_idxevent: modify according the ratio of fs_new = 500 and fs_old 
%       
%       lfpdata: only remain the channels in M1, SMA, PMC, thalamus (VA, VLo, VPLo) areas
%                bipolar for STN and GP channels
%                channel order change along areas
%
%       T_chnsarea: add a new depth variable
%                   channel order change along areas               

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
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');


fs_new = 500;


%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'preprod';

%% start here

files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
close all;
f = waitbar(0, ['Preprocessing lfp data...']);
for i = 1 : nfiles
    % wait bar
    waitbar(i/nfiles,f,['Preprocessing lfp data in file ' num2str(i) '/' num2str(nfiles)]);
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(i).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    disp(filename)

    % extract the date of the exp and the cond
    idx = strfind(filename, '_bktdt');
    dateofexp = datetime(filename(idx-6: idx-1), 'InputFormat', 'MMddyy');
    cond = parsePDCondition(datenum(dateofexp), 'Bug');
    
    % add the depth information to T_chnsarea
    T_chanDepth = T_chanDepth_extract(cond);
    T_chnsarea = addDepths2chnsarea(T_chnsarea, T_chanDepth, dateofexp);
    
    if isempty(T_chnsarea)
        disp('ignore this file as can not add new depth variable');
        
        clear lfpdata fs T_chnsarea T_idxevent
        clear filename idx dateofexp cond T_chanDepth

        continue;
    end
    
    
    % channels for m1, GM, STN and GP
    chnsM1 = find(T_chnsarea.brainarea == "M1");
    chnsPMC = find(T_chnsarea.brainarea == "PMC");
    chnsSMA = find(T_chnsarea.brainarea == "SMA");
    chnsSM = find(T_chnsarea.brainarea == "Sensory Motor");
   
    % extract M1, PMC, SMA and SM data
    lfpM1 = lfpdata(chnsM1, :, :);
    T_chnsarea_M1 = T_chnsarea(chnsM1, :);
    lfpPMC = lfpdata(chnsPMC, :, :);
    T_chnsarea_PMC = T_chnsarea(chnsPMC, :);
    lfpSMA = lfpdata(chnsSMA, :, :);
    T_chnsarea_SMA = T_chnsarea(chnsSMA, :);
    lfpSM = lfpdata(chnsSM, :, :);
    T_chnsarea_SM = T_chnsarea(chnsSM, :);
    
    
    % extract data in thalamus 
    chnsTha = find(T_chnsarea.brainarea == "VA" | T_chnsarea.brainarea == "VLo" | T_chnsarea.brainarea == "VPLo");
    lfpTha = lfpdata(chnsTha, :, :);
    T_chnsarea_Tha = T_chnsarea(chnsTha, :);
    

    % bipolar STN and GP
    chnsSTN = find(T_chnsarea.brainarea == "STN");
    lfpSTN = lfpdata(chnsSTN,:,:);
    lfpSTN = diff(lfpSTN, [],1);
    T_chnsarea_STN = T_chnsarea(chnsSTN, :);
    T_chnsarea_STN(end,:) = [];
    for chi = 1: height(T_chnsarea_STN)
        T_chnsarea_STN(chi,:).notes{1} = ['chn' num2str(chi + 1) ' - chn' num2str(chi)];
    end
    
    
    chnsGP = find(T_chnsarea.brainarea == "GP");
    lfpGP = lfpdata(chnsGP,:,:);
    lfpGP = diff(lfpGP, [],1);
    T_chnsarea_GP = T_chnsarea(chnsGP, :);
    T_chnsarea_GP(end,:) = [];
    for chi = 1: height(T_chnsarea_GP)
        T_chnsarea_GP(chi,:).notes{1} = ['chn' num2str(chi + 1) ' - chn' num2str(chi)];
    end

    
    
    
    % combine m1, thalamus and dbs channels
    lfpdata = cat(1, lfpM1, lfpPMC, lfpSMA, lfpSM, lfpTha, lfpSTN, lfpGP);
    

    T_chnsarea = cat(1, T_chnsarea_M1, T_chnsarea_PMC, T_chnsarea_SMA, T_chnsarea_SM, T_chnsarea_Tha, T_chnsarea_STN, T_chnsarea_GP);
    % change the chni 
    T_chnsarea.chni = [1:height(T_chnsarea)]';
    
    
    
    % down sample
    fs_old = fs;
    lfpdata = resampleSeg(lfpdata, fs_old, fs_new);
    T_idxevent{:,:} = T_idxevent{:,:} /fs_old * fs_new;
    fs = fs_new;
        
    
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    
    clear lfpdata fs T_chnsarea T_idxevent 
    clear idx tmpn savefilename  
    clear filename dateofexp cond T_chanDepth
    clear nM1 nSTN nGP nGM
    clear lfpm1 T_chnsarea_m1 lfpGM T_chnsarea_GM lfpSTN T_chnsarea_STN lfpGP T_chnsarea_GP
end
close(f)
disp("Processing all files Done!")
