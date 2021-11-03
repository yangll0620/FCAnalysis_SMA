function m4_seg2ShortSegments()
%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'shortSeg';


%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm3_restData_rmChns_avgArea');

% time length 
twin = 0.2;

%% Code Start Here
cond_cell = cond_cell_extract(animal);
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
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
    savefilename = [animal '_' savefilename_addstr '_' pdcond '.mat'];
    save(fullfile(savefolder, savefilename), 'lfpdata', 'fs', 'T_chnsarea');
    
    clear pdcond files
    clear lfpdata fs fs_unit T_chnsarea 
    clear savefilename
end

