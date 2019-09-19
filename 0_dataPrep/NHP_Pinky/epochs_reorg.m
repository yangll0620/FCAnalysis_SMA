function epochs_reorg()
animal = 'Pinky';
downsample_extchns(animal)

function downsample_extchns(animal)
% 1. down sample files with sample rate 3.0518e+3 to 1.073e+3
% 2. extract 1-96 utah array, 101-132 gray matter and 1-14 dbs channels

if nargin < 1
    animal = 'Pinky';
end

load('unifs.mat')
unifs = unique(fstable.fs); % unifs[1]: 1.0173e+3, unifs[2]: 3.0518e+3
n1 = 4;
n12 = round(unifs(2)/unifs(1));
if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('F:', 'yang7003@umn');
end
folder = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data', animal);
loadfolder = fullfile(folder,'epochs');
savefolder = fullfile(folder, 'epochs_reorg');
files = extractfield(dir(fullfile(loadfolder,'*.mat')), 'name');

chantbl_varNames = {'chni', 'electype', 'area', 'notes'};
chantbl_size = [96+32+14, length(chantbl_varNames)];
chantbl_varTypes = {'uint8','string','string', 'string'};
chninf = uint8([1:96+32+14]);
electype = cell(96+32+14,1);
area = cell(96+32+14,1);
notes = cell(96+32+14,1);
electype(1:96,1) = {'Utah Array'};
electype(97:128,1) = {'Gray Matter'};
electype(129:142,1) = {'DBS'};
area(1:96,1) = {'M1'};
area(97:128,1) = {'Thalamus &SMA'};
notes(1:96,1) = {''};
notes(97:128,1) = {''};

for filei = 3: length(files)
    file = files{filei};
    load(fullfile(loadfolder, file), 'lfptrial_dbs', 'lfptrial_cortical', ...
        'idxtbl_event', 'fs', 'chantbl_dbs', 'chantbl_cortical');
    
    % downsample
    if fs == unifs(1)
        n = n1;
    else if fs == unifs(2)
            n = n1 * n12;
        end
    end
    % lfptrial_cortical case
    [chn, ~, trialn] = size(lfptrial_cortical);
    for chni = 1 : chn
        tmp = squeeze(lfptrial_cortical(chni,:,:));
        tmp_downsample = downsample(tmp, n);
        if chni == 1
            tempn = size(tmp_downsample, 1);
            lfptrial = zeros(chn,tempn,trialn);
            clear tempn
        end
        lfptrial(chni, :,:) = tmp_downsample;
        clear tmp tmp_downsample
    end
    lfptrial_cortical = lfptrial;
    clear chn trialn chni lfptrial
    
    % lfptrial_dbs case
    [chn, ~, trialn] = size(lfptrial_dbs);
    for chni = 1 : chn
        tmp = squeeze(lfptrial_dbs(chni,:,:));
        tmp_downsample = downsample(tmp, n);
        if chni == 1
            tempn = size(tmp_downsample, 1);
            lfptrial = zeros(chn,tempn,trialn);
            clear tempn
        end
        lfptrial(chni, :,:) = tmp_downsample;
        clear tmp tmp_downsample
    end
    lfptrial_dbs = lfptrial;
    clear chn trialn chni lfptrial
    
    fs_new = fs /n;
    varNames = idxtbl_event.Properties.VariableNames;
    idxtbl_event = array2table(round(idxtbl_event{:,:} / fs * fs_new), 'VariableNames', varNames);
    fs = fs_new;
    clear fs_new varNames
    
        
    % extract channs
    chn = size(lfptrial_cortical, 1);
    lfptrial = cat(1,lfptrial_cortical([1:96 chn-31:chn], :,:),lfptrial_dbs);
    clear lfptrial_cortical lfptrial_dbs 
    % deal with the channel table
    area(129:142,1) = chantbl_dbs.area;
    notes(129:142,1) = chantbl_dbs.elecchn;
    
    chantbl = table('Size',chantbl_size,'VariableTypes',chantbl_varTypes,'VariableNames',chantbl_varNames);
    chantbl.chn = chninf';
    chantbl.electype = electype;
    chantbl.area = area;
    chantbl.notes = notes;
    clear chantbl_dbs chantbl_cortical
    
    % save
    save(fullfile(savefolder, file), 'lfptrial', 'fs', 'idxtbl_event', 'chantbl')
end