function [optFreezeTypes, combFreeTypes] = optFreezeTypes_extract(varargin)
%   Inputs:
%       animal
%
%       Name-Value: 
%           'codesavefolder' - code saved folder


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


optFreezeTypes = {'freeze during init Move', 'freeze during React-Reach', 'freeze during Reach', 'freeze during Manipulation'};
combFreeTypes = {'initFreeze', 'reachFreeze', 'maniFreeze'}; % combined {'freeze during React-Reach'}  and  {'freeze during Reach'} 
end