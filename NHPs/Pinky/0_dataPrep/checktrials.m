function checktrials()
% check the trials after performing SKBtrialextract(animal)
checkchnsNum()
checkfs()

function checkchnsNum()
% check the chn num of lfptrial_dbs and lfptrial_cortical 
%
animal = 'Pinky';
if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('F:', 'yang7003@umn');
end
savedir = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data');
savefolder = fullfile(savedir, animal, 'epochs');
files = extractfield(dir(savefolder), 'name');

for filei = 3: length(files)
    file = files{filei};
    load(fullfile(savefolder, file), 'lfptrial_dbs', 'lfptrial_cortical')
    
    corchns = size(lfptrial_cortical, 1);
    dbstchns = size(lfptrial_dbs, 1);
    if corchns ~= 132 || dbstchns  ~= 14
        disp([file ', cortical num = ' num2str(corchns) ', dbs num = ' num2str(dbstchns)])
    end
    clear file lfptrial_dbs lfptrial_cortical corchns dbschns
end

function  fstable = checkfs()
% check the sampling rate of each file
%
animal = 'Pinky';
if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('F:', 'yang7003@umn');
end
savedir = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data');
savefolder = fullfile(savedir, animal, 'epochs');
files = extractfield(dir(savefolder), 'name');

for filei = 3: length(files)
    file = files{filei};
    load(fullfile(savefolder, file), 'fs');
    if filei == 3
        fstable = table({file}, fs, 'VariableNames',{'file','fs'});
    else
        fstable = [fstable;table({file}, fs, 'VariableNames',{'file','fs'})];
    end
end
save('unifs.mat','fstable')

