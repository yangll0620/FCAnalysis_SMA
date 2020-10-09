function m2_STKData_narrowfiltered9_11()
%% narrow filtered STK data recorded  in frequency [9 11]Hz
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


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
tmp = char(regexp(codefilepath, 'NHPs/[A-Za-z]*', 'match'));
animal = tmp(length('NHPs/')+1:end);


%%  input setup
% band pass frequency
frebp = [9 11];

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_STKData_seg');



%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['narrowFiltered' num2str(frebp(1)) '_' num2str(frebp(2)) '-'];

%% code start here
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Narrow Filtering....']);
for filei = 1 : nfiles
    % wait bar
    waitbar(filei/nfiles,f,['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'fs', 'chnAreas');
    
    
    %%% band pass filter
    [ntemps, nchns, ntrials] = size(lfpdata);
    filterdlfp = zeros(ntemps, nchns, ntrials);
    for chni = 1 : nchns
        for triali = 1: ntrials
            filterdlfp(:,chni, triali) = filter_bpbutter(lfpdata(:,chni,triali),frebp,fs);
        end
    end
    
    
    
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ...
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'chnAreas');
    
end