function [lfpdata, fs, T_chnsarea]= seg2ShortSegments(files, twin)
if isempty(files)
    lfpdata = [];
    fs = [];
    T_chnsarea = [];
    
    return;
end

lfpdata = [];
for fi = 1: length(files)
    loadfilename = files(fi).name;
    load(fullfile(files(fi).folder, loadfilename), 'data_segments', 'fs', 'T_chnsarea');
    
    nwin = round(twin * fs);
    for segi = 1 : length(data_segments)
        seglfp = data_segments(segi).lfp;
        
        len = size(seglfp, 2);
        shortSegn = floor(len / nwin);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nwin + 1;
            endi = shortsegi * nwin;
            lfpdata = cat(3, lfpdata, seglfp(:, stri : endi));
            clear stri endi
        end
        clear seglfp len shortSegn shortsegi
    end
    
    if ~exist('fs_unit', 'var')
        fs_unit = fs;
    else
        if(fs_unit ~=fs)
            dis(['fs_unit ~=fs for file ' loadfilename])
        end
    end
    
    if ~exist('T_chnsarea_unit', 'var')
        T_chnsarea_unit = T_chnsarea;
    else
        if(~isequal(T_chnsarea_unit, T_chnsarea))
            dis(['T_chnsarea_unit ~= T_chnsarea for file ' loadfilename])
        end
    end
    
    clear nwin segi
    clear('data_segments', 'fs', 'T_chnsarea')
    clear loadfilename
    
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;


