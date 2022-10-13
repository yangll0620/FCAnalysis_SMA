function copyfile2aifolder(fromfolder, tofolder, fileformat)

if ~exist(tofolder, 'dir')
    mkdir(tofolder)
end


files =  dir(fullfile(fromfolder, ['*.' fileformat]));
for fi = 1 : length(files)
    file = files(fi).name;
    copyfile(fullfile(fromfolder, file), fullfile(tofolder, file));
    clear file
end

