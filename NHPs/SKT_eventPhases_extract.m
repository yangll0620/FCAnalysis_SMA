function EventPhases = SKT_eventPhases_extract(animal, varargin)

p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

if strcmpi(animal, 'Kitty')
    EventPhases = {'preMove'; 'earlyReach'; 'peakV'; 'lateReach'};
end

if strcmpi(animal, 'Jo')
    EventPhases = {'preMove'; 'earlyReach';  'lateReach'};
end



