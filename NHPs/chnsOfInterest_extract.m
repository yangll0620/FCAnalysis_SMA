function chnsOfI = chnsOfInterest_extract(animal, varargin)
%   
%   Usage:
%       chnsOfI = chnsOfInterest_extract('Kitty', 'codesavefolder', codesavefolder)
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



chnsOfI = {};
if strcmpi(animal,'jo')
    chnsOfI = {'M1', 'stn1-2', 'gp3-4'};
end
if strcmpi(animal,'kitty')
    chnsOfI = {'M1', 'stn1-2', 'gp1-2'};
end
end