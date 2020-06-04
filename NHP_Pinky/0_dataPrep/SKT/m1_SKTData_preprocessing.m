function m1_SKTData_preprocessing()
%% preprocessing skt data 
%
%   Processing steps as follows:
%   
%       only keep the good channels in m1 that are same as thosed used in rest(chans_m1.mat)
%
%       remove the channels not used in Gray Matter
%
%       bipolar for STN and GP channels



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
% load chans_m1.mat for channels that are used in m1
load('chans_m1.mat')


% input folder: extracted raw STK data 
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');


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
    


    % nchannels for m1, GM, STN and GP
    nM1 = length(find(T_chnsarea.brainarea == "M1"));
    nSTN = length(find(T_chnsarea.brainarea == "STN"));
    nGP = length(find(T_chnsarea.brainarea == "GP"));
    nGM = height(T_chnsarea) - nM1 - nSTN - nGP;


    % only remain data of channels in chans_m1
    lfpm1 = lfpdata(chans_m1, :, :);
    T_chnsarea_m1 = T_chnsarea(chans_m1, :);

    % remove data of band channels in GrayMatter
    T_chnsarea_GM = T_chnsarea(nM1 + 1: nM1 + nGM, :);
    lfpGM = lfpdata(nM1 + 1: nM1 + nGM, :,:);
    lfpGM = lfpGM(T_chnsarea_GM.brainarea ~= "",:,:);
    T_chnsarea_GM = T_chnsarea_GM(T_chnsarea_GM.brainarea ~= "",:);

    % bipolar STN and GP
    lfpSTN = lfpdata(nM1 + nGM + 1: nM1 + nGM + nSTN,:,:);
    lfpSTN = diff(lfpSTN, [],1);
    T_chnsarea_STN = T_chnsarea(nM1 + nGM + 1: nM1 + nGM + nSTN -1, :);
    for chi = 1: height(T_chnsarea_STN)
        T_chnsarea_STN(chi,:).notes{1} = ['chn' num2str(chi + 1) ' - chn' num2str(chi)];
    end

    lfpGP = lfpdata(nM1 + nGM + nSTN + 1: end,:,:);
    lfpGP = diff(lfpGP, [],1);
    T_chnsarea_GP = T_chnsarea(nM1 + nGM + nSTN + 1: end-1, :);
    for chi = 1: height(T_chnsarea_GP)
        T_chnsarea_GP(chi,:).notes{1} = ['GP chn' num2str(chi + 1) ' - chn' num2str(chi)];
    end

    % combine m1, GM and dbs channels
    lfpdata = cat(1, lfpm1, lfpGM, lfpSTN, lfpGP);
    T_chnsarea = cat(1, T_chnsarea_m1, T_chnsarea_GM, T_chnsarea_STN, T_chnsarea_GP);
    % change the chni 
    T_chnsarea.chni = [1:height(T_chnsarea)]';
    
    
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    
    clear lfpdata fs T_chnsarea T_idxevent 
    clear idx tmpn savefilename  
    clear filename
    clear nM1 nSTN nGP nGM
    clear lfpm1 T_chnsarea_m1 lfpGM T_chnsarea_GM lfpSTN T_chnsarea_STN lfpGP T_chnsarea_GP
end

