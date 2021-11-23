function noisy_chns = noisy_chns_extract(animal, varargin)
%
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
    noisy_chns = {};
end
if strcmpi(animal,'kitty')
    noisy_chns = {'gp3-4', 'gp4-5'};
end
if strcmpi(animal,'pinky')
    noisy_chns = {};
end
if strcmpi(animal,'bug')
    noisy_chns = {'stn0-1', 'stn4-5', 'gp5-6'};
end
end