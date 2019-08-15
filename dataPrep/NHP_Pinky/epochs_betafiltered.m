function epochs_betafiltered(animal)
% beta band filtered data, three beta bands [13 16], [16 20]Hz, and [20 30]Hz
%
% output:
%      NMRC_umn/Projects/FCAnalysis/data/Pinky /epochs_reorg_filtered/betaBand


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
savefolder = fullfile(animalfolder, 'epochs_reorg_filtered','betaBand');
files = extractfield(dir(loadfolder),'name');

% beta1, beta2, beta3 and beta
beta1 = [13 16];
beta2 = [16 20];
beta3 = [20 30];
beta = cat(1,beta1, beta2, beta3);

filessaved = extractfield(dir(savefolder),'name');

% skip the first 2, '.', '..'
for i = 3: length(files)
    file = files{i};
    
    disp(['dealing i = ' num2str(i) ':' file])
    
    if ~ismember(file, filessaved)
        
        % load epochs file
        load(fullfile(loadfolder, file), 'lfptrial', 'fs', 'idxeventtbl','chantbl', 'idxevent' , 'idxevent_varNames');
        
        % filtered in beta1
        [lfptrial_beta1]= betafiltered(lfptrial, beta1, fs);
        
        % filtered in beta2
        [lfptrial_beta2]= betafiltered(lfptrial, beta2, fs);
        
        % filtered in beta3
        [lfptrial_beta3]= betafiltered(lfptrial, beta3, fs);
        
        % save the filtered data
        save(fullfile(savefolder, file), 'lfptrial_beta1','lfptrial_beta2','lfptrial_beta3', 'beta', ...
            'fs', 'idxeventtbl','chantbl','idxevent', 'idxevent_varNames');
    end
end


function data_betafiltered = betafiltered(data, betaband, fs)
% data: chn * tempn * trialn

[chn, ~, trialn] = size(data);
data_betafiltered = zeros(size(data));
for chni = 1: chn
    for triali = 1: trialn
        x = squeeze(data(chni, :, triali));
        data_betafiltered(chni, :, triali)= filter_bpbutter(x,betaband,fs);
    end
end

    