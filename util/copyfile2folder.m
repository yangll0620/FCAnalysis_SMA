function copyfile2folder(codefilepath, savefolder)

[~, codefilename]= fileparts(codefilepath);

status = copyfile([codefilepath '.m'], fullfile(savefolder, [codefilename '.m']));
disp(['copied code status = ' num2str(status)])