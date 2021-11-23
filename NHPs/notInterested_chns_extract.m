function notAOI_chns = notInterested_chns_extract(animal, varargin)
%
%   Input:
%       aimal:
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
    notAOI_chns = {};
end
if strcmpi(animal,'kitty')
    notAOI_chns = {};
end
if strcmpi(animal,'pinky')
    notAOI_chns = {'lCd', 'lSMA', 'lVA', 'lVLo', 'lVPLo', 'rMC', 'rSMA', 'rVA', 'rVLo', 'rVPLo'};
end
if strcmpi(animal,'bug')
    notAOI_chns = {};
end
end