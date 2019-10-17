function codecorresfolder = code_corresfolder(codefilepath, makefolder, makesubfolder)
%% code_corresfolder 
%           return the corresponding folder for codefilepath
%       and create the corresponding folder if not exist (makefolder = True)
%		and create the subfolders 'out' and 'store' subfolders
%
%   Usage:
%           codecorresfolder = code_corresfolder(codefilepath, true)
% 
%   args:
% 
%       codefilepath: the full path and the name of code file without suffix 
%                     (i.e codefilepath = '/home/lingling/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/util/folder_extract for folder_extract.m)


if nargin < 3
    makesubfolder = false;
end


if nargin < 2
    makefolder = true;
end

[codefolder, codefilename] = fileparts(codefilepath);

% extract the substring of codefolder after 'code'
subfolder = codefolder(strfind(codefolder, 'code') + length('code'):end);

% delete the first character if the first one is '/'
if subfolder(1) == '/'
    subfolder(1) = [];
end

% extract pipeline folder
[~, ~, pipelinefolder, ~] = exp_subfolders();

% return the corresponding folder in pipeline folder
codecorresfolder = fullfile(pipelinefolder, subfolder, codefilename);

% make folder if needed
if makefolder == true && ~exist(codecorresfolder, 'dir')
    mkdir(codecorresfolder)
end


% make sub folders if needed
if makesubfolder == true 

	% out subfolder
	outsubfolder = fullfile(codecorresfolder, 'out');
	if ~exist(outsubfolder, 'dir')
    	mkdir(outsubfolder)
    end


   	% store subfolder
	storesubfolder = fullfile(codecorresfolder, 'store')
	if ~exist(storesubfolder, 'dir')
    	mkdir(storesubfolder);
    end

end