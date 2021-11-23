function [clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal)
if strcmpi(animal,'jo')
%     clim_Spectrogram_STN = [-130 -110];
%     clim_Spectrogram_GP = [-130 -110];
%     clim_Spectrogram_Others = [-120 -90];
    clim_Spectrogram_STN = [-30 -15];
    clim_Spectrogram_GP = [-30 -15];
    clim_Spectrogram_Others = [-30 -15];
end
if strcmpi(animal,'kitty')
    %clim_Spectrogram_STN = [-130 -100];
    %clim_Spectrogram_GP = [-130 -100];
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