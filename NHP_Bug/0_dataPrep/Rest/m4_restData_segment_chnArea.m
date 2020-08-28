function m4_restData_segment_chnArea()
%% seg into intervals with same time length (remove bad segments) and add chn-area inf
%
%   Processing steps as follows:
%       1. segment into intervals with same time length (remove bad segments)
%
%       2. concatenate lfpseg_GM, lfpseg_stn and lfpseg_gp into lfp
%       
%       3. add chn-area information T_chnsarea using file  'Bug_GMChnAreaInf.csv'


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx1 = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx1 + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% folder
[datafolder, ~, ~, ~] = exp_subfolders();
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);



%%  input setup
segt = 2;

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_M1Power');


% gray matter chn-area information file
filename_GMChnsarea = [animal '_GMChnAreaInf.csv'];
file_GMChnsarea =  fullfile(datafolder, animal, filename_GMChnsarea);
nGM = 96;nSTN = 7; nGP = 7;


%% save setup
savefolder = codecorresfolder;
savefilename_prefix = [animal 'seged_'];
saveformat_datestr = 'yyyymmdd';

%% Code Start Here
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Seging Rest Data and Adding Chn-Area Inf....']);


for filei = 1 : nfiles
    disp(files(filei).name)
    
    % wait bar
    waitbar(filei/nfiles,f,['Seging Rest Data and adding Chn-Area Inf in file ' num2str(filei) '/' num2str(nfiles)]);
      
    
    filename = files(filei).name;
    filefolder = files(filei).folder;
    
    
    % date of exp, bktdt
    idx1 = strfind(filename, '_tdt');
    dateofexp = datenum(filename(idx1-8:idx1-1), 'yyyymmdd');
    idx2 = strfind(filename, '.mat');
    tdtbk = str2num(filename(idx1+4:idx2-1));
    
    pdcondition = parsePDCondition(dateofexp, animal);
    
    
    %%% seg data  %%%
    [lfpsegs_GM, lfpsegs_stn, lfpsegs_gp, fs]  = segPerFile(fullfile(filefolder, filename), segt);
    
    
    %%% concatenate  %%%
    lfpdata = cat(2, lfpsegs_GM, lfpsegs_stn, lfpsegs_gp);
    
    
    %%% add chn-area information T_chnsarea  %%%
    [T_chnsarea_GM] = chanInf_GM(file_GMChnsarea, nGM, pdcondition);
    [T_chnsarea] = chanInf_addDBStoGM(T_chnsarea_GM, nSTN, nGP);
    
    
    
    % save
    savefilename = [savefilename_prefix  pdcondition '_' datestr(dateofexp, saveformat_datestr) '_bktdt' num2str(tdtbk)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea');
    
    
    
    
    clear filefolder filename pdcondition
    clear idx1 dateofexp idx2 tdtbk
    clear lfpsegs_GM lfpsegs_stn lfpsegs_gp fs
    clear lfpdata T_chnsarea_GM T_chnsarea
end
end



function [lfpsegs_GM, lfpsegs_stn, lfpsegs_gp, fs]  = segPerFile(file_data, segt)
%
% args:
%   file_data: the file containing data_segments, segRemain and fs (full path, e.g ../m2_restData_selectSeg_M1Power/Bug_cleanedRestData_mild_20190206_tdt1.mat)
%   segt: the seg durantion (s)
% 
% returns:
%   lfpsegs_GM, lfpsegs_stn, lfpsegs_gp: the segmented lfp data (ntemp * nchns * ntrials)
%   fs: sample rate


load(file_data, 'data_segments', 'segsRemain', 'fs');

lfpsegs_GM =[];
lfpsegs_stn =[];
lfpsegs_gp =[];
for segi = 1: size(data_segments,2)
    
    if segsRemain(segi) ==0 % bad segment, not used
        continue;
    end
    
    seg_lfpGM = data_segments(segi).lfp_array;
    seg_lfpstn = data_segments(segi).lfp_stn;
    seg_lfpgp = data_segments(segi).lfp_gp;
    
    
    % seg
    ntemp = size(seg_lfpGM,1);
    nsegtemp = ceil(fs * segt);
    tsegn = floor(ntemp/nsegtemp);
    for i = 1: tsegn
        idx_str = (i-1) * nsegtemp + 1;
        idx_end = i * nsegtemp;
        
        lfpsegs_GM = cat(3, lfpsegs_GM, seg_lfpGM(idx_str:idx_end,:));
        lfpsegs_stn = cat(3, lfpsegs_stn, seg_lfpstn(idx_str:idx_end,:));
        lfpsegs_gp = cat(3, lfpsegs_gp, seg_lfpgp(idx_str:idx_end,:));
        
        clear idx_str idx_end
    end
    
    clear seg_lfpm1 seg_lfpGM seg_lfpstn seg_lfpgp
    clear ntemp nsegtemp tsegn
end
end


function [T_chnsarea] = chanInf_GM(file_GMChnsarea, nGM, pdcondition)
% extract gray matter channel inf table
%   
%   Args:
%       file_GMCchnsarea: the file storing the gray matter chn-area inf
%       nGM: the total number of gray matter channels
%       pdcondition: pd condition
%
%   Output:
%       T_chnsarea: table of Gray matter channel inf,
%                   (T_chnsarea.Properties.VariableNames: 
%                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})


% initial the attitudes of chanarea table T_chnsarea: 
% chni_vec, electype, brainareas, notes, recordingchn
chni_vec = uint8([1: nGM]');
electype = cell(nGM,1);
brainareas = cell(nGM,1);
notes = cell(nGM,1);
recordingchn = zeros(nGM,1);
electype(1:nGM,1) = {'Gray Matter'}; % electype



% deal with Gray Matter
brainareas(1:nGM,1) = {''};
T = readtable(file_GMChnsarea);

if strcmp(pdcondition, 'normal')
    T.channels = T.channels_normal;
end

if strcmp(pdcondition, 'mild')
    T.channels = T.channels_mild;
end

for i = 1 : length(T.brainarea)
    area = T.brainarea{i};
    tmpcell = split(T.channels{i}, ',');
    
    for j = 1 : length(tmpcell)
        chn = str2num(char(tmpcell{j}));
        brainareas(chn,1) = {area};
        
        clear chn
    end
end
recordingchn(1:nGM) = [1:nGM]';
notes(1:nGM,1) = {''};


% channel information table
T_chnsarea = table;
T_chnsarea.chni = chni_vec;
T_chnsarea.brainarea = brainareas;
T_chnsarea.recordingchn = recordingchn;
T_chnsarea.electype = electype;
T_chnsarea.notes = notes;

end


function [T_chnsarea] = chanInf_addDBStoGM(T_chnsarea_GM, nSTN, nGP)
% add DBS channel to the end of the table GM inf
%   
%   Args:
%       T_chnarea_GM: the chn-area inf table for gray matter
%       nstn, ngp: the number of stn and gp channels
%
%   Output:
%       T_chnsarea: table of Gray matter + DBS channel inf,
%                   (T_chnsarea.Properties.VariableNames: 
%                   {'chni'}  {'brainarea'}    {'recordingchn'}    {'electype'}    {'notes'})


nGM = height(T_chnsarea_GM);


T_chnsarea = T_chnsarea_GM;
for chi = 1: nSTN
    T_chnsarea = [T_chnsarea; {nGM + chi, 'STN', chi, 'DBS', ''}];
end
for chi = 1: nGP
    T_chnsarea = [T_chnsarea; {nGM + nSTN + chi, 'GP', chi, 'DBS', ''}];
end


end
