function m2_STKData_Downsample()
%% downsample STK Data
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

% the corresponding pipeline and the parent folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup
fs_new = 500;

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_preprocessing');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['downsampled'];



%% Start Here

files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['dowmsampling Data....']);

for filei = 1 : nfiles
    disp(files(filei).name)
    
    % wait bar
    waitbar(filei/nfiles,f,['downsample STK Data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data
    filename = files(filei).name;
    load(fullfile(inputfolder, filename),  'fs', 'lfpdata', 'T_chnsarea','T_idxevent', 'chans_m1', 'GMChnAreas');
    
    
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
    save(fullfile(savefolder, savefilename), 'lfpdata','fs','T_chnsarea','T_idxevent', 'chans_m1', 'GMChnAreas');
    
    
end
end

function resampledX = resampleSeg(X, fs_old, fs_new)
%% downsample along the second dim X:  nchns * ntemp * ntrials
%
%   return:
%       resampledX: nchns * ntemp(resampled) * ntrials

resampledX = [];
for triali = 1: size(X, 3)
    x = X(:,:,triali); % x: nchns * ntemp 
    resamp = resample(x', round(fs_new), round(fs_old)); % resamp: ntemp(resampled) * nchns
    resamp = resamp'; % resamp: nchns * ntemp(resampled)
    
    resampledX = cat(3, resampledX, resamp); % resampledX: nchns * ntemp(resampled) * ntrials
end

end