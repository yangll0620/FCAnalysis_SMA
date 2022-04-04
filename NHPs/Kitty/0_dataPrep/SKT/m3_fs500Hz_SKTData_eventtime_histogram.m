function m3_fs500Hz_SKTData_eventtime_histogram()
%  extract lfp data respect to reachonset
% 
%  return:
%        lfptrials: nchns * ntemp * ntrials


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

%% global variables

% animal
animal = animal_extract(codecorresfolder);


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

%%  input setup
inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI_goodReach');

pdcond = 'moderate';

% saved fig format
savefig_format = 'tif';

%% Code Start Here

% load T_idxevent from all files
files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
T_idxevent = table();
for fi = 1: length(files)
    load(fullfile(inputfolder, files(fi).name), 'fs_ma', 'T_idxevent_ma');
    if ~exist('fs','var')
        fs = fs_ma;
    else
        if fs_ma ~= fs
            disp(['fs ~= fs_ma for ' files(fi).name]);
            continue;
        end
    end
    T_idxevent =[ T_idxevent; T_idxevent_ma];
    clear fs_ma T_idxevent_ma 
end
clear files

% calc reaction time, manipulating time
reactime = (T_idxevent.ReachTimeix - T_idxevent.TargetTimeix) / fs;
reachime = (T_idxevent.TouchTimeix - T_idxevent.ReachTimeix) / fs;
maniputime = (T_idxevent.ReturnTimeix - T_idxevent.TouchTimeix) / fs;


% hist plot
nbins = 30;
subplot(3, 1, 1);
histogram(reactime, nbins);
title(['reaction time, ntrials = ' num2str(length(reactime))])
subplot(3, 1, 2);
histogram(reachime, nbins);
title(['reach time, ntrials = ' num2str(length(reachime))])
subplot(3, 1, 3);
histogram(maniputime, nbins);
title(['manipulate time, ntrials = ' num2str(length(maniputime))])

savefile = fullfile(savefolder, ['eventtime_hist.' savefig_format]);
saveas(gcf, savefile);

