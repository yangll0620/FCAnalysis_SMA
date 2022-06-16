function freqsOfI = freqsOfInterest_extract(animal, varargin)
%   frequency bands interested 
%   Usage:
%       freqsOfI = chnsOfInterest_extract('Kitty', 'codesavefolder', codesavefolder)
%   
%   Inputs:
%       animal
%
%       Name-Value: 
%           'codesavefolder' - code saved folder


% parse params
p = inputParser;
addParameter(p, 'ishighfreq', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

ishighfreq = p.Results.ishighfreq;
codesavefolder = p.Results.codesavefolder;

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


switch animal
    case 'Kitty'
        fhigh_AOIs = [220 240];
    case 'Jo'
        fhigh_AOIs = [230 250
                      320 340;];
    otherwise
        fhigh_AOIs = [];
end

if ishighfreq
    freqsOfI = fhigh_AOIs;
else
    freqsOfI = [8 40];
end

end