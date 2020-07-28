function m3_STKData_narrowfiltered19_21()
%% narrow filtered STK data recorded  in frequency [29 31]Hz
%
%   based on spectrogram and psd analysis in
%   pipeline/../m2_SKTData_usefulChannel


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup
% band pass frequency
frebp = [19 21];
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_usefulChannel');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['narrowFiltered' num2str(frebp(1)) '_' num2str(frebp(2)) '-'];

%% starting: narrow filter the lfp data of all the files
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Narrow Filtering....']);
tic;
for filei = 1 : nfiles
    elapsedTime = toc;
    
    % wait bar
    waitbar(filei/nfiles,f,['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles) ...
         ', elapseTime ' num2str(elapsedTime)]);
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    
    % band pass filter
    [nchns, ntemps, ntrials] = size(lfpdata);
    filterdlfp = zeros(nchns, ntemps, ntrials);
    for chni = 1 : nchns
        for triali = 1: ntrials
            filterdlfp(chni,:, triali) = filter_bpbutter(lfpdata(chni,:,triali),frebp,fs);
        end
    end
         
    % convert T_chnsarea, T_idxevent into Cell/Matrix
    [chnsarea_Cell, chnsarea_vNames]= conv_tbl2Cell(T_chnsarea);
    [idxevent_Matrix, idxevent_vNames]= conv_tbl2matrix(T_idxevent);
    
    % save
    lfpdata = filterdlfp;
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'chnsarea_Cell', 'chnsarea_vNames', 'idxevent_Matrix', 'idxevent_vNames');
    
    
    clear lfpdata fs T_chnsarea
    clear idx tmpn savefilename preproc_lfp
    clear areas filename
    
end
close(f);
disp(['narrow filtered lfpdata  are saved to ' savefolder])

