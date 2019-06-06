function epochs_betafiltered(animal)
if nargin < 1
    animal = 'Pinky';
end

addpath(fullfile('..', '..', 'util'))

if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('F:', 'yang7003@umn');
end
folder = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data', animal);
loadfolder = fullfile(folder,'epochs_reorg');
savefolder = fullfile(folder, 'epochs_reorg_filtered');
files = extractfield(dir(loa
dfolder),'name');
beta1 = [13 16];
beta2 = [16 20];
beta3 = [20 30];
beta = cat(1,beta1, beta2, beta3);
filessaved = extractfield(dir(savefolder),'name');
for i = 3: length(files)
    file = files{i};
    if ~ismember(file, filessaved)
        load(fullfile(loadfolder, file), 'lfptrial', 'fs', 'idxtbl_event','chantbl');
        [lfptrial_beta1]= betafiltered(lfptrial, beta1, fs);
        [lfptrial_beta2]= betafiltered(lfptrial, beta2, fs);
        [lfptrial_beta3]= betafiltered(lfptrial, beta3, fs);
        save(fullfile(savefolder, file), 'lfptrial_beta1','lfptrial_beta2','lfptrial_beta3','fs', 'idxtbl_event','chantbl','beta');
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

    