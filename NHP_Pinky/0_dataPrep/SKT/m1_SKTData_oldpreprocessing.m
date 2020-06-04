function m1_SKTData_oldpreprocessing()
%% preprocessing skt data 
%
%   Processing steps as follows:
%       down sampled to 500Hz
%   
%       common average reference inside one area (except DBS) deleted
%       
%       band pass filtered in frebp = [2 100];


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup
% band pass frequency
frebp = [2 100];
% downsampled to fs_new Hz
fs_down = 500;
% input folder: extracted raw rest data 
inputfolder = fullfile(codecorresParentfolder, 'm0_SKTData_extract');


%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'preprod';


%% Start:  preprocess the lfp data of all the files
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
close all;
f = waitbar(0, ['Preprocessing lfp data...']);
for i = 1 : nfiles
    % wait bar
    waitbar(i/nfiles,f,['Preprocessing lfp data in file ' num2str(i) '/' num2str(nfiles)]);
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(i).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea', 'T_idxevent');
    
    % get all the areas
    areas = {};
    for j = 1 : height(T_chnsarea)
        area = T_chnsarea.brainarea{j};
        if sum(cellfun(@(x) strcmp(x,  area), areas)) == 0
            areas = [areas, area];
        end
        clear area
    end
    
    
    % preprocessing lfp data in each area
    for j = 1: length(areas)
        area = areas{j};
        chns_area = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
        
        % lfp: nchns * ntemp * ntrials
        lfp = lfpdata(chns_area, :, :);
        
        % preprocess lfp
        [tmp] = downsample_lfp(lfp, fs, fs_down);
%         if ~(strcmp(area, 'STN') || strcmp(area, 'GP'))
%             % do CARRef for the areas not DBS leads
%             tmp = CARRef_lfp(tmp);
%         end
        tmp = filtered_lfp(tmp, frebp, fs_down);
        
        if j == 1
            [nchns, ~, ntrials] = size(lfpdata);
            preproc_lfp = zeros(nchns, size(tmp,2), ntrials);
        end
        preproc_lfp(chns_area, :, :) = tmp;
         
        clear area chns_area tmp
    end
    
    % save
    lfpdata = preproc_lfp;
    T_idxevent = adjust_Idxevent(T_idxevent, fs, fs_down);
    fs = fs_down;    
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea', 'T_idxevent');
    
    
    clear lfpdata fs T_chnsarea T_idxevent 
    clear idx tmpn savefilename preproc_lfp T_idxevent_adjust
    clear areas filename
end
close(f)


function [donwsampled_lfp] = downsample_lfp(lfp, fs, fs_down)
%% downsample lfp data 
% 
%   Input
%       lfp: lfp data [nchns ntemp ntrials]
%       fs: sample rate of the lfp data
%       fs_down: new sample rate that is downsampled to
%       T_idxevent: idxevent table, will be changed once downsampled
%       
%   Return
%       donwsampled_lfp: downsampled lfp data [nchns ntemp_new ntrials]
%       T_idxevent: the adjusted T_idxevent accroding to the down sample

[nchns, ~, ntrials] = size(lfp);

for triali = 1: ntrials
    tmp = squeeze(lfp(:, :,triali)); % tmp: nchns * ntemps
    
    % down sample, treats each column as a separate sequence.
    tmp = tmp'; % tmp: ntemps * nchns
    downsampled_tmp = resample(tmp,round(fs_down),round(fs)); % downsampled_tmp :  ntemps * nchns
    downsampled_tmp = downsampled_tmp'; % downsampled_tmp :   nchns * ntemps
    
    
    % assignment
    if triali == 1
        donwsampled_lfp = zeros(nchns, size(downsampled_tmp,2), ntrials);
    end
    donwsampled_lfp(:, :,triali) = downsampled_tmp;
    
end

function T_idxevent_adjust = adjust_Idxevent(T_idxevent, fs, fs_down)
% adjust the T_idxevent
varNames = T_idxevent.Properties.VariableNames;
T_idxevent_adjust = array2table(round(T_idxevent{:,:} / fs * fs_down), 'VariableNames', varNames);



function refed_lfp = CARRef_lfp(lfp)
%% common average reference lfp data 
%
%   Input
%       lfp: lfp data [nchns ntemp ntrials]
%       
%   Return
%       refed_lfp: common average referenced lfp data [nchns ntemp ntrials]


% common average reference in one area
refed_lfp = lfp - repmat(mean(lfp, 1), size(lfp,1),1,1);


function filterd_lfp = filtered_lfp(lfp, frebp, fs)
%% band pass filter lfp data 
% 
%   Input
%       lfp: lfp data [nchns ntemp ntrials]
%       frebp: the band pass frequency
%       fs: sample rate of the lfp data
%       
%   Return
%      filterd_lfp: filetered lfp data [nchns ntemp ntrials]


[nchns, ntemp, ntrials] = size(lfp);
filterd_lfp = zeros(nchns, ntemp, ntrials);

% band pass filter
for chni = 1 : nchns
    for triali = 1: ntrials
        tmp = squeeze(lfp(chni, :, triali));
        filterd_lfp(chni, :, triali) = filter_bpbutter(tmp,frebp,fs);
        clear tmp
    end
end
    
