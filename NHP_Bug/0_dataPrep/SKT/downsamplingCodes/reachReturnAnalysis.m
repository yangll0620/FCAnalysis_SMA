function [reachTimes_normal, returnTimes_normal, reachTimes_mild, returnTimes_mild]= reachReturnAnalysis()

%% input
inputfolder = '.';


cond = 'normal';
files = dir(fullfile(inputfolder, ['*_' cond '_*.mat']));
[reachTimes, returnTimes]= reachReturnTimes_fromAllfiles(files);
reachTimes_normal = reachTimes;
returnTimes_normal = returnTimes;
clear reachTimes returnTimes

cond = 'mild';
files = dir(fullfile(inputfolder, ['*_' cond '_*.mat']));
[reachTimes, returnTimes]= reachReturnTimes_fromAllfiles(files);
reachTimes_mild = reachTimes;
returnTimes_mild = returnTimes;
clear reachTimes returnTimes



function [reachTimes, returnTimes]= reachReturnTimes_fromAllfiles(files)
% extract all the reach and return times from all the files
%
%   Return:
%       reachTimes, returnTimes: ntrials (the trial number from all the files) * 1
%   

reachTimes = [];
returnTimes = [];
for filei = 1: length(files)
    file = fullfile(files(filei).folder, files(filei).name);
    
    load(file, 'T_idxevent', 'fs');
    
    
    reachTime = (T_idxevent.TouchTimeix - T_idxevent.ReachTimeix)/fs;
    returnTime = (T_idxevent.MouthTimeix - T_idxevent.ReturnTimeix)/fs;
    
    
    reachTimes = [reachTimes; reachTime];
    returnTimes = [returnTimes; returnTime];
    
    clear file T_idxevent fs
    clear reachTime returnTime
end


