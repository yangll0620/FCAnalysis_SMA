function SKBtrialextract(animal)
% SKBtrialextract extract all the single Kluver board task for animal
%  down sampled to 500Hz
% saveto: 
%      ['H:\My Drive\NMRC_umn\Projects\FunctionalConnectivityAnalysis\data\' animal]
% 
% Description:
%     deal with the PD condition ('normal', 'mild', or 'moderate')

%% input parameters
if nargin < 1
    animal = 'Pinky';
end

googlesheetlink_Pinky = '1Mn_HcvWt4FVc2kcvRMbDj5jQffyGNrwppNOK6T5W-zI';

if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('H:', 'My Drive');
end
savedir = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data');
savefolder = fullfile(savedir, animal);
if ~exist(savefolder)
    mkdir(savefolder)
end

%% variable name for each column in google sheet
varname_taskdisp = 'Brief Description';
varname_tdtbk = 'TDT Block';
varname_date = 'OutputFolderName';
taskdisps = {'Kluver', 'Single Target Kluver', 'SKB', 'Single KB'  ...
    'Single Kluver', 'Single', 'Single ' ,'Single Target', ...
    'Single Target Kluver', 'Single-Target Kluver', 'Single-target Kluver',...
    'single', 'single Kluver', 'single-target Kluver'};% different discriptions of Kluver board task

%% add util path
addpath(genpath(fullfile('..','..', 'util')))

%% get the google sheet
googlesheetdata = GetGoogleSpreadsheet(googlesheetlink_Pinky);
varname_googlesheet = googlesheetdata(1,:);
col_match = @(varname) find(cell2mat(cellfun(@(x) contains(x, varname),varname_googlesheet,'UniformOutput', 0)));
coli_taskdisp = col_match(varname_taskdisp);
coli_tdtbk = col_match(varname_tdtbk);
coli_date =  col_match(varname_date);

strs_taskdispinsheet = googlesheetdata(1:end,coli_taskdisp);
row_match = @(strpattern) find(cell2mat(cellfun(@(x) ~isempty(find(strcmp(strpattern, x))), strs_taskdispinsheet, 'UniformOutput',0)));
idxrow = row_match(taskdisps); % the index for one particular task
if ispc
    NMRCdriver = 'Y:';
end
if isunix
    NMRCdriver = '/run/user/1000/gvfs/ftp:host=nmrc_dserver1.local/root2';
end
      
processedfolder = fullfile(NMRCdriver, 'Animals2', 'Pinky', 'Recording', 'Processed', 'DataDatabase');

% check whether date of that day has been processed
for i = 1:length(idxrow)
    rowi = idxrow(i);
    block = googlesheetdata{rowi, coli_tdtbk};
    datefoldername = googlesheetdata{rowi, coli_date}; % datename : 'Pinky_012017'
    tmp = split(datefoldername, '_');
    dateofexp = datenum(tmp{2},'mmddyy'); % 
    onedaypath = fullfile(processedfolder, datefoldername);
    
    % ma file
    mafileexisttag = 0;
    mafolder = fullfile(processedfolder, datefoldername, ['Block-' num2str(block)]);
    mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
    if ~isempty(mafilestruct)
        mafileexisttag = 1;
    end
    
    % LFP file original saved in .sevfile
    lfpexisttag = 0;
    folder_lfp = fullfile(onedaypath, 'LFP', ['Block-' num2str(block)]);
    files = dir(folder_lfp);
    if ~isempty(files)
        filename = extractfield(files, 'name');
        match = cellfun(@(x) ~isempty(regexp(x, ['LFPch[0-9]*.nex'],'match')),filename,'UniformOutput',false);
        if ~isempty(find(cell2mat(match))) % exist *LFPchn#.nex file
            lfpexisttag = 1;
        end
    end
    
    % DBSLFP
    lfpdbsexisttag = 0;
    dbslfpfolder = fullfile(onedaypath, 'DBSLFP', ['Block-' num2str(block)]); % Z:\Animals\Jo\Recording\Processed\DataDatabase\Pinky_101315\DBSLFP\Block-8
    dbslfpfile = fullfile(dbslfpfolder, [animal '_GrayMatter_eyetracking_DT1_' ...
        datestr(dateofexp, 'mmddyy') '_Block-' num2str(block) '_DBSLFP.nex']);
    if ~isempty(dbslfpfile)
        lfpdbsexisttag = 1;
    end
    
    if mafileexisttag && lfpexisttag && lfpdbsexisttag
        
        disp(['The ' num2str(i) 'th, rowi = ' num2str(rowi) ': ' animal '_' datestr(dateofexp, 'mmddyy') '-Block' num2str(block)])
        
        [lfptrial_cortical, lfptrial_dbs,fs,idxtbl_event ,chantbl_cortical, chantbl_dbs] = extractlfptrial(onedaypath, block);
        
        % save
        condition = parsePDCondtion_Pinky(dateofexp);
        savefile = [animal '_lfptrial_' condition '_' datestr(dateofexp,'mmddyy') '_block' num2str(block) '.mat'];
        save(fullfile(savefolder, savefile), 'lfptrial_cortical', 'lfptrial_dbs','fs','idxtbl_event','chantbl_cortical', 'chantbl_dbs')

    end
end
end

function lfptrial_new = resample_lfptrial(lfptrial, fs_now, fs_new)
% Input:
%   lfptrial: chns * ntemporal * ntrials
[chns, ~, ntrials] = size(lfptrial);
for chni = 1: chns
    for triali = 1: ntrials
        x = squeeze(lfptrial(chni,:, triali));
        y = resample(x,fs_new,fs_now);
        
        if chni == 1 && triali == 1
            ntemporal_new = length(y);
            lfptrial_new = zeros(nchns, ntemporal_new ,ntrial);
        end
        lfptrial_new(chni, :, triali) =y;
        clear tmp y
    end
end
end
