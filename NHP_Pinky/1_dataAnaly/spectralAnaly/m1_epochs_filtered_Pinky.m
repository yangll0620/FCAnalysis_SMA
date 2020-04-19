function epochs_filtered_Pinky_1(freq)
%  band filtered data
%   
% Args:
%   freq: the frequency band
%
% output:
%      corresponding folder/out


%% codecorresfolder
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));


[~, ~, pipelinefolder, ~] = exp_subfolders();

codecorresfolder = code_corresfolder(codefilepath, true, true);


%% load and save folders
loadfolder = fullfile(pipelinefolder, '0_dataPrep', 'epochs');
savefolder = fullfile(codecorresfolder, 'out');


filessaved = extractfield(dir(fullfile(savefolder)),'name');

% filter each file
files = extractfield(dir(fullfile(loadfolder, '*.mat')),'name');
for i = 1: length(files)
    file = files{i};
    
    disp(['dealing i = ' num2str(i) ': ' file])
    
    if ~ismember(file, filessaved)
        
        % load epochs file
        load(fullfile(loadfolder, file), 'lfptrial', 'fs', 'idxeventtbl','chantbl', 'idxevent' , 'idxevent_varNames');
        
        % filtered in fre
        [lfptrial_fre]= frefiltered(lfptrial, freq, fs);
        
        % save the filtered data
        savefile = ['freq[' num2str(freq(1)) ' ' num2str(freq(2)) ']' file]; 
        save(fullfile(savefolder, savefile), 'lfptrial_fre', 'freq', ...
            'fs', 'idxeventtbl','chantbl','idxevent', 'idxevent_varNames');
        
    end
end


function data_frefiltered = frefiltered(data, freband, fs)
% data: chn * tempn * trialn

[chn, ~, trialn] = size(data);
data_frefiltered = zeros(size(data));
for chni = 1: chn
    for triali = 1: trialn
        x = squeeze(data(chni, :, triali));
        data_frefiltered(chni, :, triali)= filter_bpbutter(x,freband,fs);
    end
end

    