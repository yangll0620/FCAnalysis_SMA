function m2_restData_preprocessing()
%% preprocessing rest data recorded
%
%   Processing steps as follows:
%       down sampled to 500Hz
%   
%       common average reference inside one area
%       
%       band pass filtered in frebp = [1 100];


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
% band pass frequency
frebp = [1 100];
% downsampled to fs_new Hz
fs_new = 500;
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm1_restData_extract');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'preprod';


%% Start:  preprocess the lfp data of all the files
files = dir(fullfile(inputfolder, '*.mat'));
n = length(files);
f = waitbar(0, ['Preprocessing']);
for i = 1 : n
    % wait bar
    waitbar(i/n,f,['Preprocessing lfp data in file ' num2str(i) '/' num2str(n)]);
    
    filename = files(i).name;
    
    % lfpdata: nchns * ntemp
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea');
    
    % transpose lfpdata
    lfpdata = lfpdata';
    
    % get all the areas
    areas = {};
    for j = 1 : height(T_chnsarea)
        area = T_chnsarea.brainarea{j};
        if sum(cellfun(@(x) strcmp(x,  area), areas)) == 0
            areas = [areas, area];
        end
        clear area
    end
    
    % preprocessing lfp data in each area
    for j = 1: length(areas)
        area = areas{j};
        chns_area = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
        lfp = lfpdata(:, chns_area);
        
        % preprocess lfp
        tmp = preprocessing_lfp(lfp, frebp, fs, fs_new);
        if j == 1
            [~, nchns] = size(lfpdata);
            preproc_lfp = zeros(size(tmp,1),nchns);
        end
        preproc_lfp(:,chns_area) = tmp;
        
        clear area chns tmp
    end
    fs = fs_new;   
    
    % save
    lfpdata = preproc_lfp;
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea');
    
    
    clear lfpdata fs T_chnsarea
    clear idx tmpn savefilename preproc_lfp
    clear areas filename
    
end
close(f);
disp(['preprocessed lfpdata using grayMatter are saved to ' savefolder])

function preproc_lfp = preprocessing_lfp(lfp, frebp, fs, fs_new)
%% preprocessing lfp data (downsample + common average reference + band pass filtered)
% 
%   Input
%       lfp: lfp data [ntemp nchns]
%       frebp: the band pass frequency
%       fs: sample rate of the lfp data
%       fs_new: new sample rate that is downsampled to
%       
%   Return
%       


% down sample, treats each column as a separate sequence.
lfp = resample(lfp,round(fs_new),round(fs));


% common average reference in one area
reflfp = lfp - repmat(mean(lfp, 2), 1, size(lfp,2));

% band pass filter
filterdreflfp = zeros(size(reflfp));
for chi = 1 : size(reflfp, 2)
    filterdreflfp(:, chi) = filter_bpbutter(reflfp(:,chi),frebp,fs_new);
end

preproc_lfp = filterdreflfp;
