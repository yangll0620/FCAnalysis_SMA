function STKtrialextract()
% STKtrialextract extract all the single Kluver board task marked by Ying
% for Pinky
%  
% 1.  read the stk information in
% /Projects/FCAnalysis/metainf/pinky_skbinf.csv and only extract the
% dateset marked 'Yes' in column YingUsed

% 2. Inputs:
%       preprocessed MA data in root2:\Animals2\Pinky\Recording\Processed\DataDatabase

% 3.  down sample files with sample rate 3.0518e+3 to 1.073e+3
%
% 4. extract 1-96 utah array from LFP\Block-1\*_LFPchn1-chn96.nex and 101-132 gray matter LFP\Block-1\*_LFPchn101-chn132.nex and 1-14 dbs
% channels from DBSLFP\Block-1\*_DBSLFP.nex
%
% 5. STK(stk) is short for Single Target Kluver
%
% saveto: 
%      ['\NMRC_umn\Projects\FunctionalConnectivityAnalysis\data\' animal]
% 
% Description:
%     deal with the PD condition ('normal', 'mild', or 'moderate')

%% animal related information
animal = 'Pinky';

% channel number of M1, SMA, DBS
nM1 = 96; nSMA = 32; nDBS = 14;

%% add util path
addpath(genpath(fullfile('..','..','util')))

%% various folders extraction
% set google drive and root2 for unix and windows separately
if isunix
    NMRC_umn = fullfile('/home', 'lingling', 'yang7003@umn.edu', 'NMRC_umn');
    root2 = '/home/lingling/root2';
else
    if ispc
        NMRC_umn = fullfile('F:', 'yang7003@umn', 'NMRC_umn');
        root2 = 'Y:';
    end
end

% folder for the processed data in root2
folder_root2processed = fullfile(root2, 'Animals2', animal, 'Recording', 'Processed', 'DataDatabase');

% folder for FCAnalysis project /NMRC_umn/Projects/FCAnalysis
folder_FCProject = fullfile(NMRC_umn, 'Projects','FCAnalysis');

%% Read stk information in /Projects/FCAnalysis/metainf/pinky_skbinf.csv
% metainf folder for animal
folder_metainf = fullfile(folder_FCProject, 'metainf', animal);
% metainf file name
inffilename = 'pinky_skbinf.csv';
inffile = fullfile(folder_metainf, inffilename);

% open the text file
fid = fopen(inffile, 'r');

% extract the column name
varNames = split(fgetl(fid), ',');

% Format for each line of text:
%   column1: test (%s)
%	column2: int8 (%d8)
%   column3: int8 (%d8)
%	column4: categorical (%C)
%   column5: categorical (%C)
%	column6: categorical (%C)
formatSpec = '%s%d8%d8%s%s%s%[^\n\r]';
delimiter = ',';
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndofLine', '\r\n');

fclose(fid);

% Create output variable
pinkyskbinf = table(dataArray{1:end-1}, 'VariableNames', varNames);

% Clear temporary variables
clearvars folder_metainf inffilename inffile delimiter startRow formatSpec fid dataArray varNames;


%% extract all STK trials
% folder for segmented epochs
folder_epochsdata = fullfile(folder_FCProject, 'data', animal,'epochs');
% create folder for segmented epochs if not exist
if exist(folder_epochsdata, 'file') ~= 7
    mkdir(folder_epochsdata)
end

% extract the dateofexp, bkma, and bktdt of used STK trials which are marked
% 'Yes' in the column 'YingUsed'
validskbs = pinkyskbinf{strcmp(pinkyskbinf.YingUsed, 'Yes'), {'dateofexp','bkma', 'bktdt'}};

% extract all STK trials
for i = 1 : length(validskbs)
    
    % date of experiment
    dateofexp = datenum(validskbs(i, 1), 'yymmdd');
    % tdt block number
    bktdt = char(validskbs(i, 3));
    
    
    % one day path
    onedaypath = fullfile(folder_root2processed, [animal '_' datestr(dateofexp, 'mmddyy')]);
    
    disp(['extracting ' onedaypath '-tdtbk' num2str(bktdt)])
    
    % extract trials of lfp data for particular date and tdt block number
    [lfptrial_cortical, lfptrial_dbs,fs,idxeventtbl, chantbl_dbs] = extractlfptrial(onedaypath, bktdt);
    
    if isempty(lfptrial_cortical) || isempty(lfptrial_dbs) || isempty(fs) || isempty(idxeventtbl) ||isempty(chantbl_dbs)
        continue;
    end
        
    % concatenate along the first channel dim, get only the 96 chns M1
    % and 32 chns SMA
    chn_cortical = size(lfptrial_cortical, 1);
    lfptrial = cat(1,lfptrial_cortical([1:96 chn_cortical-31:chn_cortical], :,:),lfptrial_dbs);

    % down sample to around 250 Hz
    [lfptrial_down, fs_down] = lfp_downsample(lfptrial, fs);
    lfptrial = lfptrial_down;
    
    % get the event time inf and channel inf
    [idxeventtbl, chantbl] = chaninf(idxeventtbl, chantbl_dbs, nM1, nSMA, nDBS, fs, fs_down);
    
    % convert table idxtbl_event animalto matrix and varNames as python cann't
    % read matlab table structure
    [idxevent, idxevent_varNames]= conv_tbl2matrix(idxeventtbl);
    
    % set the sample rate 
    fs = fs_down;
    
    % get the pd conditioon for the date of experiment
    pdcondition = parsePDCondtion_Pinky(dateofexp);
    
    % name of the saved file 
    savefile = [animal '_lfptrial_' pdcondition '_' datestr(dateofexp,'mmddyy') '_bktdt' num2str(bktdt) '.mat'];
    
    % save
    save(fullfile(folder_epochsdata, savefile), 'lfptrial','fs','idxeventtbl','chantbl','idxevent_varNames','idxevent')
    
    
    clear dateofexp bktdt onedaypath pdcondition savefile; 
    clear lfptrial_cortical lfptrial_dbs fs idxtbl_event chantbl_dbs  chn_cortical lfptrial;
    clear lfptrial_down fs_down idxtbl_event chantbl idxevent idxevent_varNames;  
end
end


function [lfptrial_down, fs_down] = lfp_downsample(lfptrial, fs)
%
% down sample files with sample rate 1.073e+3, 3.0518e+3 to 254.31Hz (1.073e+3 / 4)
% 
% Inputs:
%       lfptrial: nchns * ntimes * ntrials
% 
% Output:
%       lfptrial_down: downsampled lfptrial (nchns * ntimes *ntrials)
%       fs_down: the new down sampling rate
%       

load('unifs.mat','unifs');
n1 = 4;
n12 = round(unifs(2)/unifs(1));

% downsample
if fs == unifs(1)
    n = n1;
else
    if fs == unifs(2)
        n = n1 * n12;
    end
end

% 
[chn, ~, trialn] = size(lfptrial);
for chni = 1 : chn
    tmp = squeeze(lfptrial(chni,:,:));
    tmp_downsample = downsample(tmp, n);
    if chni == 1
        tempn = size(tmp_downsample, 1);
        lfptrial_down = zeros(chn,tempn,trialn);
        clear tempn
    end
    lfptrial_down(chni, :,:) = tmp_downsample;
    clear tmp tmp_downsample
end
fs_down = fs/n;
end

function [idxeventtbl, chantbl] = chaninf(idxeventtbl, chantbl_dbs, nM1, nSMA, nDBS, fs, fs_new)
% extract the event idx inf table and channel inf table
%
%   Output:
%       idxeventtbl: table of event idx inf
%       chantbl: table of channel inf
%


chninf = uint8([1:nM1+nSMA+nDBS]);
electype = cell(nM1+nSMA+nDBS,1);
area = cell(nM1+nSMA+nDBS,1);
notes = cell(nM1+nSMA+nDBS,1);
electype(1:nM1,1) = {'Utah Array'};
electype(nM1+1:nM1+nSMA,1) = {'Gray Matter'};
electype(nM1+nSMA+1:nM1+nSMA+nDBS,1) = {'DBS'};
area(1:nM1,1) = {'M1'};
area(nM1+1:nM1+nSMA,1) = {'Thalamus &SMA'};
notes(1:nM1,1) = {''};
notes(nM1:nM1+nSMA,1) = {''};

varNames = idxeventtbl.Properties.VariableNames;
idxeventtbl = array2table(round(idxeventtbl{:,:} / fs * fs_new), 'VariableNames', varNames);

% deal with the DBS channel table
area(nM1+nSMA+1:nM1+nSMA+nDBS,1) = chantbl_dbs.area;
notes(nM1+nSMA+1:nM1+nSMA+nDBS,1) = chantbl_dbs.elecchn;

% channel information table
chantbl = table;
chantbl.chn = chninf';
chantbl.electype = electype;
chantbl.area = area;
chantbl.notes = notes;
end

function [matrix, varNames]= conv_tbl2matrix(tbl)
% python could not read the matlab table objects, thus convert the
% table structure to matrix and varNames, table variables should be in the
% same type

matrix = tbl{:,:};
varNames = tbl.Properties.VariableNames;

end
