function status = copyfile2folder(codefilepath, savefolder)

if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

[~, codefilename]= fileparts(codefilepath);

status = copyfile([codefilepath '.m'], fullfile(savefolder, [codefilename '.m']));