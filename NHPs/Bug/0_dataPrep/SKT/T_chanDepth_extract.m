function T_chanDepth = T_chanDepth_extract(cond)
% extract the table T_chanDepth 
% 
%   Arg:
%      cond: the pd condition ('mild' or 'normal')
%
%  Return:
%       T_chanDepth: table for channel depth with brain area as varnames
%       E.g.
%           T_chanDepth(1:3,1:5)
% 
%               3×5 table
% 
%                   Var1       M1      M1_1     M1_2      M1_3 
%                 ________    _____    ____    ______    ______
% 
%                      NaT       63     55         73        70
%                 02/27/19    26.25     50     33.125    20.625
%                 03/05/19    26.25     50     33.125    20.625



%% find datafolder
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

[datafolder, ~, ~, ~] = exp_subfolders();

%% Start Here
if strcmp(cond, 'mild')
    filename_chanDepth = 'Bug_channelDepth_mild.csv';
else
    if strcmp(cond, 'normal')
        filename_chanDepth = 'Bug_channelDepth_normal.csv';
    end
end    
file_chanDepth = fullfile(datafolder, filename_chanDepth);

% return the import options 
opts = detectImportOptions(file_chanDepth);

nvars = length(opts.VariableNames);

% read format for each column
fmt = ['%{MM/dd/yy}D' repmat('%f', 1, nvars-1)];

% read the channel depth information
warning off
T_chanDepth = readtable(file_chanDepth, 'Format',fmt);
warning on