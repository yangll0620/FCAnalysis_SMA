function [clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal, varargin)
%   Inputs:
%       animal
%
%       Name-Value: 
%           'F_AOI': frequency of Interested, default [8 40]
%
%           'codesavefolder' - code saved folder

p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'f_AOI', [8 40], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
f_AOI = p.Results.f_AOI;

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

clim_Spectrogram_STN = [];
clim_Spectrogram_GP = [];
clim_Spectrogram_Others =[];

if isequal(f_AOI, [8 40])
    if strcmpi(animal,'jo')
        clim_Spectrogram_STN = [-30 -15];
        clim_Spectrogram_GP = [-30 -15];
        clim_Spectrogram_Others = [-30 -10];
    end
    if strcmpi(animal,'kitty')
        clim_Spectrogram_STN = [-40 0];
        clim_Spectrogram_GP = [-40 0];
        clim_Spectrogram_Others = [-40 0];
    end
    if strcmpi(animal,'bug')
        clim_Spectrogram_STN = [-130 -90];
        clim_Spectrogram_GP = [-140 -120];
        clim_Spectrogram_Others = [-130 -90];
    end
    if strcmpi(animal,'pinky')
        clim_Spectrogram_STN = [-35 -15];
        clim_Spectrogram_GP = [-35 -15];
        clim_Spectrogram_Others = [-35 -15];
    end
end

if isequal(f_AOI, [40 400])
    if strcmpi(animal,'jo')
        clim_Spectrogram_STN = [-45 -25];
        clim_Spectrogram_GP = [-45 -25];
        clim_Spectrogram_Others = [-45 -25];
    end
    if strcmpi(animal,'kitty')
        clim_Spectrogram_STN = [-50 -20];
        clim_Spectrogram_GP = [-50 -20];
        clim_Spectrogram_Others = [-50 -20];
    end
end


end