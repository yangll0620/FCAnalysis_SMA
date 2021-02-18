function m2_STKData_narrowfiltered9_11()
%% 
%
% input:
%   m1_SKTData_avgArea
%
% output:
%   skip the file with only one trial


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
[fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
animal = codecorresfolder(fi + length('NHPs/'):j);


%%  input setup
% band pass frequency
frebp = [9 11];
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_avgArea');

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
    waitbar(filei/nfiles,f,[' [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles) ...
         ', elapseTime ' num2str(elapsedTime)]);
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    if(height(T_idxevent) < 10)
        disp([filename ' has less than 10 trials, skip!']);
        continue;
    end
    
    
    % band pass filter
    [nchns, ntemps, ntrials] = size(lfpdata);
    filterdlfp = zeros(nchns, ntemps, ntrials);
    
    for chni = 1 : nchns
        for triali = 1: ntrials
            filterdlfp(chni,:, triali) = filter_bpbutter(lfpdata(chni,:,triali),frebp,fs);
        end
    end
         
    % extract chnAreas cell for used in python
    chnAreas =T_chnsarea.brainarea;
    idx_STN = find(strcmp(T_chnsarea.brainarea, 'STN'));
    for i = 1: length(idx_STN)
        chnAreas{idx_STN(i)} = ['stn' num2str(i-1) '-' num2str(i)];
    end
    idx_GP = find(strcmp(T_chnsarea.brainarea, 'GP'));
    for i = 1: length(idx_GP)
        chnAreas{idx_GP(i)} = ['gp' num2str(i-1) '-' num2str(i)];
    end

    % extract idxevent Matrix for used in python
    idxevent = T_idxevent{:,:};

    
    % save
    lfpdata = filterdlfp;
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'chnAreas', 'idxevent', 'T_chnsarea', 'T_idxevent');
    
    
    clear lfpdata fs T_chnsarea
    clear idx tmpn savefilename preproc_lfp
    clear areas filename
    
end
close(f);
disp(['narrow filtered lfpdata  are saved to ' savefolder])

