function checktrials()

animal = 'Pinky';
if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('H:', 'My Drive');
end
savedir = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data');
savefolder = fullfile(savedir, animal);
files = extractfield(dir(savefolder), 'name');

for filei = 3: length(files)
    file = files{filei};
    load(fullfile(savefolder, file), 'lfptrial_dbs', 'lfptrial_cortical')
    
    corchns = size(lfptrial_cortical, 1);
    dbstchns = size(lfptrial_dbs, 1);
    if corchns ~= 132 || dbstchns  ~= 14
        disp([file 'cortical num = ' num2str(corchns) ', dbs num = ' num2str(dbstchns)])
    end
    clear file lfptrial_dbs lfptrial_cortical corchns dbschns
end
    
    

