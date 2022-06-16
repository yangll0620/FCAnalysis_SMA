function m2_restData_selectSeg_Power_fs1000Hz()
%   marked the good and bad segments
%
%   Processing steps as follows:
%       1. mark the segment good (1) using fs500 data
%
%
%       2. bipolar for STN and GP channels
%
%       3. Down sample trials into fs_new = 1000
%
%       4. combine lfp_m1, lfp_stn and lfp_gp into lfp(ntemp * nchns)
%
%       5. add table T_chnsarea (nchns * 1  cell)


%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code') - 1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder, 'util')));
addpath(genpath(fullfile(codefolder, 'NHPs')));

% the corresponding pipeline folder for this code
[codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

%% global parameters

%% input setup
folder_fs500 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Jo\0_dataPrep\Rest';

infolder_fs500_befsel = fullfile(folder_fs500, 'm1_restData_cleaned_extract');
infolder_fs500_aftsel = fullfile(folder_fs500, 'm2_restData_selectSeg_Power');

fs_new = 1000;

%% save setup
savefolder = codecorresfolder;
codesavefolder = fullfile(savefolder, 'code');

%% Code Start Here

% copy code into
copyfile2folder(mfilename('fullpath'), codesavefolder);

files = dir(fullfile(infolder_fs500_aftsel, '*.mat'));
nfiles = length(files);
for fi = 1:nfiles  
    filename = files(fi).name;

    file_befsel = fullfile(infolder_fs500_befsel, filename);
    file_fs500_aftsel = fullfile(infolder_fs500_aftsel, filename);

    % load
    load(file_befsel, 'fs', 'data_segments');
    data_aftsel_fs500 = load(file_fs500_aftsel);


    if length(data_segments) == length(data_aftsel_fs500.data_segments) %
        selectSegs = ones(length(data_segments),1);
    else
        continue;
    end

    for segi = 1:length(data_segments)
        % lfp_m1: ntemp * nchns
        lfp_m1 = data_segments(segi).lfp_m1;

        % bipolar lfp_stn and lfp_gp
        lfp_stn = diff(data_segments(segi).lfp_stn, [], 2);
        lfp_gp = diff(data_segments(segi).lfp_gp, [], 2);

        %%% add chn-area information T_chnsarea  %%%
        if ~exist('T_chnsarea', 'var')
            nM1 = size(lfp_m1, 2);
            nSTN = size(lfp_stn, 2);
            nGP = size(lfp_gp, 2);
            T_chnsarea = chanInf_M1DBS(nM1, nSTN, nGP);

            clear nM1 nSTN nGP
        end


        %%% down sample %%%
        lfp_m1 = resample(lfp_m1, round(fs_new), round(fs));
        lfp_stn = resample(lfp_stn, round(fs_new), round(fs));
        lfp_gp = resample(lfp_gp, round(fs_new), round(fs));

        %%% combine lfp_m1, lfp_stn and lfp_gp %%%%
        data_segments(segi).lfp = cat(2, lfp_m1, lfp_stn, lfp_gp);
        clear lfp_m1 lfp_stn lfp_gp
    end

    fs = fs_new;

    % remove lfp_m1, lfp_stn and lfp_gp as they have been combined together
    data_segments = rmfield(data_segments, {'lfp_m1', 'lfp_stn', 'lfp_gp'});

    savefile = fullfile(savefolder, filename);
    save(savefile, 'fs', 'data_segments', 'T_chnsarea', 'selectSegs');

    clear filename  file_befsel  file_fs500_aftsel
    clear data_aftsel_fs500 fs data_segments T_chnsarea selectSegs
    
end

end

function [T_chnsarea] = chanInf_M1DBS(nM1, nSTN, nGP)
% add M1 and DBS channel together
%
%   Args:
%       nM1: the number of M1 channels
%       nstn, ngp: the number of stn and gp channels
%
%   Output:
%       T_chnsarea: table of M1 + DBS channel inf,
%                   (T_chnsarea.Properties.VariableNames:
%                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})

% initial the attitudes of chanarea table T_chnsarea:
% chni_vec, electype, brainareas, notes, recordingchn
chni_vec = uint8([1:nM1]');
electype = cell(nM1, 1);
brainareas = cell(nM1, 1);
notes = cell(nM1, 1);
recordingchn = uint8([1:nM1]');
electype(1:nM1, 1) = {'Gray Matter'}; % electype
brainareas(1:nM1, 1) = {'M1'}; % electype

% channel information table of M1
T_chnsarea = table;
T_chnsarea.chni = chni_vec;
T_chnsarea.brainarea = brainareas;
T_chnsarea.recordingchn = recordingchn;
T_chnsarea.electype = electype;
T_chnsarea.notes = notes;

% add DBS lead
for chi = 1:nSTN
    T_chnsarea = [T_chnsarea; {nM1 + chi, 'STN', chi, 'DBS', ''}];
end

for chi = 1:nGP
    T_chnsarea = [T_chnsarea; {nM1 + nSTN + chi, 'GP', chi, 'DBS', ''}];
end

end
