function status = copyfile2folder(codefilepath, savefolder)

[~, codefilename]= fileparts(codefilepath);

status = copyfile([codefilepath '.m'], fullfile(savefolder, [codefilename '.m']));
if status
    disp(['copied ' codefilename ' .m to ' savefolder])
end