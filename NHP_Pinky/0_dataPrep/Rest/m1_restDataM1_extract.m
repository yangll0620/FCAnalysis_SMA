function m1_restDataM1_extract()
%% extract the rest data with same dateofexp in skb analysis in Area M1

%% codecorresfolder
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

[datafolder, ~, pipelinefolder, ~] = exp_subfolders();

codecorresfolder = code_corresfolder(codefilepath, true, false);

%% preprocessed folder in root2
processedfolder_inroot2 = '/home/lingling/root2/Animals2/Pinky/Recording/Processed/DataDatabase';


%% save setup
animal = 'Pinky';
savefolder = codecorresfolder;
savefilename_prefix = [animal '_restdataM1_'];

%% read restinf
file_restinf= fullfile(pipelinefolder, 'NHP_Pinky', '0_dataPrep', '0_restinf', 'Pinky_restinf.csv');
data_restinf = readtable(file_restinf);

chns_M1 = [1:96];
for i = 1: length(data_restinf.dateofexp)
    
    % date of exp, bktdt, pdcond for one rest
    dateofexp = datetime(num2str(data_restinf.dateofexp(i)), 'InputFormat','yyMMdd');
    bktdt = data_restinf.bktdt(i);
    pdcond = data_restinf.pdCondition{i};
    
    nexfolder_restlfp = fullfile(processedfolder_inroot2, ['Pinky_' datestr(dateofexp, 'mmddyy')],...
        'LFP', ['Block-' num2str(bktdt)]);
    
    % extract the prefix of the files and chns in folder_restlfp
    [file_prefix, ~] = nexfilenames_infolder(nexfolder_restlfp);
    
    %
    [lfpdata, fs] = lfpdataextract_fromnexfileofallchns(nexfolder_restlfp, file_prefix, chns_M1);
    
    if ~isempty(lfpdata)
        
        
        savefile = fullfile(savefolder, [savefilename_prefix  pdcond '_' ...
            datestr(dateofexp, 'mmddyy') '_bktdt' num2str(bktdt)]);
        save(savefile, 'lfpdata', 'fs');
    end
    
    
    clear dateofexp bktdt pdcond nexfolder_restlfp file_prefix chns i
end