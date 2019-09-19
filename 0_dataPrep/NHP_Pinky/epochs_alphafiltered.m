function epochs_alphafiltered(animal)
% alpha band filtered data
%
% output:
%      NMRC_umn/Projects/FCAnalysis/data/Pinky /epochs_filtered/alphaBand


if nargin < 1
    animal = 'Pinky';
end

addpath(fullfile('..', '..', 'util'))

% set google drive and root2 for unix and windows separately
if isunix
    NMRC_umn = fullfile('/home', 'lingling', 'yang7003@umn.edu', 'NMRC_umn');
else
    if ispc
        NMRC_umn = fullfile('F:', 'yang7003@umn', 'NMRC_umn');
    end
end

animalfolder = fullfile(NMRC_umn, 'Projects','FCAnalysis','data', animal);
loadfolder = fullfile(animalfolder,'epochs');
savefolder = fullfile(animalfolder, 'epochs_filtered','alphaBand');
files = extractfield(dir(loadfolder),'name');

% alpha
alpha = [8 13];

filessaved = extractfield(dir(savefolder),'name');

% skip the first 2, '.', '..'
for i = 3: length(files)
    file = files{i};
    
    disp(['dealing i = ' num2str(i) ':' file])
    
    if ~ismember(file, filessaved)
        
        % load epochs file
        load(fullfile(loadfolder, file), 'lfptrial', 'fs', 'idxeventtbl','chantbl', 'idxevent' , 'idxevent_varNames');
        
        % filtered in alpha
        [lfptrial_alpha]= alphafiltered(lfptrial, alpha, fs);
        
        % save the filtered data
        savefile = ["alphaband_" file]; 
        save(fullfile(savefolder, savefile), 'lfptrial_alpha', 'alpha', ...
            'fs', 'idxeventtbl','chantbl','idxevent', 'idxevent_varNames');
        
    end
end


function data_alphafiltered = alphafiltered(data, alphaband, fs)
% data: chn * tempn * trialn

[chn, ~, trialn] = size(data);
data_alphafiltered = zeros(size(data));
for chni = 1: chn
    for triali = 1: trialn
        x = squeeze(data(chni, :, triali));
        data_alphafiltered(chni, :, triali)= filter_bpbutter(x,alphaband,fs);
    end
end

    