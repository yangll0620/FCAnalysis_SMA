function SKBtrialextract(animal)
% SKBtrialextract extract all the single Kluver board task for animal
%  
% 1.  find the date of that day has been processed

% 2.  down sample files with sample rate 3.0518e+3 to 1.073e+3
%
% 3. extract 1-96 utah array, 101-132 gray matter and 1-14 dbs channels
%
% saveto: 
%      ['\NMRC_umn\Projects\FunctionalConnectivityAnalysis\data\' animal]
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
    googledrive = fullfile('F:', 'yang7003@umn');
end
savedir = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data');
savefolder = fullfile(savedir, animal, 'epochs');
if ~exist(savefolder)
    mkdir(savefolder)
end

%% variable name for each column in google sheet
varname_taskdisp = 'MPTP injection'; % Brief Description column record task discription
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
% varname_googlesheet: column name of google sheet
varname_googlesheet = googlesheetdata(1,:);

% extract the column number for 'Brief Description', 'tdtbk' and 'date'
col_match = @(varname) find(cell2mat(cellfun(@(x) contains(x, varname),varname_googlesheet,'UniformOutput', 0)));
coli_taskdisp = col_match(varname_taskdisp);
coli_tdtbk = col_match(varname_tdtbk);
coli_date =  col_match(varname_date);

% extract the task discription
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
        
        [lfptrial_cortical, lfptrial_dbs,fs,idxtbl_event, chantbl_dbs] = extractlfptrial(onedaypath, block);
        
        % concatenate along the first channel dim, get only the 96 chns M1
        % and 32 chns SMA
        chn_cortical = size(lfptrial_cortical, 1);
        lfptrial = cat(1,lfptrial_cortical([1:96 chn_cortical-31:chn_cortical], :,:),lfptrial_dbs);
        clear chn_cortical
        
        % down sample to around 250 Hz
        [lfptrial_down, fs_down] = lfp_downsample(lfptrial, fs);
        lfptrial = lfptrial_down;
        
        % get the event time inf and channel inf
        nM1 = 96; nSMA = 32; nDBS = 14;
        [idxtbl_event, chantbl] = chaninf(idxtbl_event, chantbl_dbs, nM1, nSMA, nDBS, fs, fs_down);
        [idxevent, idx_varNames]= conv_tbl2matrix(idxtbl_event);
        fs = fs_down;
        clear lfptrial_down fs_down
        
        % save
        condition = parsePDCondtion_Pinky(dateofexp);
        savefile = [animal '_lfptrial_' condition '_' datestr(dateofexp,'mmddyy') '_block' num2str(block) '.mat'];
        save(fullfile(savefolder, savefile), 'lfptrial','fs','idxtbl_event','chantbl','idx_varNames','idxevent')
        
        clear lfptrial_cortical lfptrial_dbs lfptrial
    end
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
else if fs == unifs(2)
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

function [idxtbl_event, chantbl] = chaninf(idxtbl_event, chantbl_dbs, nM1, nSMA, nDBS, fs, fs_new)
% extract the event idx inf table and channel inf table
%
%   Output:
%       idxtbl_event: table of event idx inf
%       chantbl: table of channel inf
%
chantbl_varNames = {'chni', 'electype', 'area', 'notes'};
chantbl_size = [nM1+nSMA+nDBS, length(chantbl_varNames)];
chantbl_varTypes = {'uint8','string','string', 'string'};

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

    
varNames = idxtbl_event.Properties.VariableNames;
idxtbl_event = array2table(round(idxtbl_event{:,:} / fs * fs_new), 'VariableNames', varNames);

% deal with the DBS channel table
area(nM1+nSMA+1:nM1+nSMA+nDBS,1) = chantbl_dbs.area;
notes(nM1+nSMA+1:nM1+nSMA+nDBS,1) = chantbl_dbs.elecchn;

chantbl = table('Size',chantbl_size,'VariableTypes',chantbl_varTypes,'VariableNames',chantbl_varNames);
chantbl.chn = chninf';
chantbl.electype = electype;
chantbl.area = area;
chantbl.notes = notes;
end

function [matrix, varNames]= conv_tbl2matrix(tbl)
% python could not read the matlab table objects, thus convert the
% table structure to matrix and varNames

matrix = tbl{:,:};
varNames = tbl.Properties.VariableNames;

end
