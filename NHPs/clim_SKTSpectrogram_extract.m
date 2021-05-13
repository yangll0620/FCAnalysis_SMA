function [clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal)
if strcmpi(animal,'jo')
    clim_Spectrogram_STN = [-130 -110];
    clim_Spectrogram_GP = [-130 -110];
    clim_Spectrogram_Others = [-120 -90];
end
if strcmpi(animal,'kitty')
    clim_Spectrogram_STN = [-130 -100];
    clim_Spectrogram_GP = [-130 -100];
    clim_Spectrogram_Others = [];
end
if strcmpi(animal,'bug')
    clim_Spectrogram_STN = [-130 -90];
    clim_Spectrogram_GP = [-130 -115];
    clim_Spectrogram_Others = [-110 -80];
end
if strcmpi(animal,'pinky')
    clim_Spectrogram_STN = [-130 -115];
    clim_Spectrogram_GP = [-130 -110];
    clim_Spectrogram_Others = [-120 -90];
end

end