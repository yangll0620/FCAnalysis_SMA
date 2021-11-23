function unwanted_DBS = unwanted_DBS_extract(animal, varargin)
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

if strcmpi(animal,'jo')
    unwanted_DBS = {'stn3-4', 'stn4-5', 'stn5-6', 'stn6-7', 'gp0-1'};
end
if strcmpi(animal,'kitty')
    unwanted_DBS = {'stn3-4', 'stn4-5', 'stn5-6', 'stn6-7', 'gp6-7'};
end
if strcmpi(animal,'pinky')
    unwanted_DBS = {'stn0-1', 'stn1-2', 'stn2-3', 'stn6-7', 'gp0-1'};
end
if strcmpi(animal,'bug')
    unwanted_DBS = {'stn5-6', 'stn6-7'};
end
end