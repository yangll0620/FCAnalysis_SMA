function [codecorresfolder, codecorresParentfolder] = output_code_corresfolder(codefilepath, varargin)
%% output_code_corresfolder 
%           code_corresfolder version for output
%
%           return the corresponding folder and the parent folder of the corresponding folder in output folder for codefilepath
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
%
%       Name-Value: 
%           'makefolder' -- making folder (true) or not (false)
%


% parse params
p = inputParser;
addParameter(p, 'makefolder', true, @(x) isscalar(x)&&islogical(x));


parse(p,varargin{:});
makefolder = p.Results.makefolder;


[codefolder, codefilename] = fileparts(codefilepath);

% extract the substring of codefolder after 'code'
subfolder = codefolder(strfind(codefolder, 'code') + length('code'):end);

% delete the first character if the first one is '/'
if subfolder(1) == '/'
    subfolder(1) = [];
end

% extract pipeline folder
[~, ~, ~, outputfolder] = exp_subfolders();

% return the corresponding folder in pipeline folder
codecorresfolder = fullfile(outputfolder, subfolder, codefilename);

% return the parent folder of the corresponding folder in pipeline folder
codecorresParentfolder = fullfile(outputfolder, subfolder);

% make folder if needed
if makefolder == true && ~exist(codecorresfolder, 'dir')
    mkdir(codecorresfolder)
end