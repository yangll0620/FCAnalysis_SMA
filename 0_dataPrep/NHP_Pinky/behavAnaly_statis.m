function behavAnaly_statis(animal)
% statistical analysis of reation, reach, manipulate and return time for
% normal and mild states

if nargin < 1
    animal = 'Pinky';
end

if isunix
    drive = fullfile('/home', 'lingling');
end

if ispc
    drive = fullfile('F:', 'yang7003@umn');
end

savedir = fullfile(drive, 'NMRC_umn','Projects','FCAnalysis','data');
savefolder = fullfile(savedir, animal, 'behavior');
if ~exist(savefolder)
    mkdir(savefolder)
end


[tdur_normal_calc, tdur_mild_calc, tdur_normal_read, tdur_mild_read, files_SKB] = behavExtract(animal);
for i = 1: 4
    [~, p_calc(i)]= ttest2(tdur_normal_calc(:,i), tdur_mild_calc(:, i));
    [~, p_read(i)]= ttest2(tdur_normal_read(:,i), tdur_mild_read(:, i));
end
avg_normal_calc = mean(tdur_normal_calc,1);
avg_mild_calc = mean(tdur_mild_calc,1);
avg_normal_read = mean(tdur_normal_read,1);
avg_mild_read = mean(tdur_mild_read,1);

varNames = {'avg_normal_calc', 'avg_mild_calc', 'p_calc' ,'avg_normal_read', 'avg_mild_read', 'p_read'};
timeNames = {'reaction_time', 'reach_time', 'manipulation_time', 'return_time'};
tbl_normild_sta = table(avg_normal_calc', avg_mild_calc', p_calc', avg_normal_read', avg_mild_read', p_read','VariableNames',varNames);
tbl_normild_sta.Properties.RowNames = timeNames;

save(fullfile(savefolder, 'normild_statisRes.mat'), 'tbl_normild_sta', 'files_SKB')
end

function [tdur_normal_calc, tdur_mild_calc, tdur_normal_read, tdur_mild_read, files_SKB] = behavExtract(animal)
% extract time duration for reaction, reach, manipulation and return subphases of the single Kluver board task for animal
%  
% load:
%       *SingleTargetKluver_Analyze2.mat in preprocessed folder
% 
% Output:
%   tdur_normal_calc: ntrials * 1
%       time duration for reaction, reach, manipulation and return
%       subphases in normal, calculated besed on idx 
%   tdur_mild_calc: ntrials * 1
%       time duration for reaction, reach, manipulation and return
%       subphases in mild, calculated besed on idx
%   tdur_normal_read: ntrials * 1
%       time duration for reaction, reach, manipulation and return
%       subphases in normal, read the precalculated time from *SingleTargetKluver_Analyze2.mat 
%   tdur_mild_read: ntrials * 1
%       time duration for reaction, reach, manipulation and return
%       subphases in mild, read the precalculated time from *SingleTargetKluver_Analyze2.mat 
%
%   files_SKB: nfiles
%       the used SKB files
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
tevent_mild = [];
tevent_normal = [];
tdur_mild_read = [];
tdur_normal_read = [];
files_SKB = ["files_SKB"];
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
        
        clear file_skb
        
        [tevent_ma, tevent_varNames] = calcMAeventTime(onedaypath, block);
        [tdur_ma, tdur_varNames] = readMAdur(onedaypath, block);
        
        [ntrials_tdur, ~] = size(tdur_ma);
        [ntrials_tevent, ~] = size(tevent_ma);
        if ntrials_tdur ~= ntrials_tevent
           disp(['The ' num2str(i) 'th, rowi = ' num2str(rowi) ': ' animal '_' datestr(dateofexp, 'mmddyy') '-Block' num2str(block)])
           size(tdur_ma)
           size(tevent_ma)
        end
        
        
        condition = parsePDCondtion_Pinky(dateofexp);

        file_skb  = [animal '_' condition '_' datestr(dateofexp, 'mmddyy') '-Block' num2str(block)];
        files_SKB = [files_SKB; file_skb];
        if strcmp(condition, 'normal')
            tevent_normal = cat(1, tevent_normal, tevent_ma);
            tdur_normal_read = cat(1, tdur_normal_read, tdur_ma);
        else
            if strcmp(condition, 'mild')
                tevent_mild = cat(1, tevent_mild, tevent_ma);
                tdur_mild_read = cat(1, tdur_mild_read, tdur_ma);
            end
        end
        
        clear tdur_ma
    end
end
[tdur_normal_calc] = tdur_calc(tevent_normal, tevent_varNames);
[tdur_mild_calc] = tdur_calc(tevent_mild, tevent_varNames);

end



function [tdur, tdur_varNames] = tdur_calc(tevent, tevent_varNames)
%
% calculate the reaction time, reach time, manipulation time and return
% time
%
% Output:
%       tdur: ntrials * 4
%
%
col_target = find(strcmp(tevent_varNames, 'target'));
col_reach = find(strcmp(tevent_varNames, 'reach'));
col_touch = find(strcmp(tevent_varNames, 'touch'));
col_return = find(strcmp(tevent_varNames, 'return'));
col_mouth = find(strcmp(tevent_varNames, 'mouth'));

cols_reaction = [col_target, col_reach];
cols_reach = [col_reach, col_touch];
cols_manipulate = [col_touch, col_return];
cols_return = [col_return, col_mouth];

tdur = [];
cols = cols_reaction;
tdur_reaction = tevent(:, cols(2)) - tevent(:,cols(1));
tdur = cat(2, tdur, tdur_reaction);
cols = cols_reach;
tdur_reach = tevent(:, cols(2)) - tevent(:,cols(1));
tdur = cat(2, tdur, tdur_reach);
cols = cols_manipulate;
tdur_manipulate = tevent(:, cols(2)) - tevent(:,cols(1));
tdur = cat(2, tdur, tdur_manipulate);
cols = cols_return;
tdur_return = tevent(:, cols(2)) - tevent(:,cols(1));
tdur = cat(2, tdur, tdur_return);
tdur_varNames = {'reaction_time', 'reach_time', 'manipulation_time', 'return_time'};
end


function [tevent_ma, tevent_varNames] = calcMAeventTime(onedaypath, block)
% calculate the MA event time slot from one file based on idx 
%
%
%  Inputs:
%   onedaypath:  'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_013017'
%   block: 10  
%
%  Outputs:
%   tevent_ma : matrix of event time slot ntrials * 5
%   tevent_varNames: names of each col in tevent_ma 1 *5 cell 
%
mafolder = fullfile(onedaypath, ['Block-' num2str(block)]); %  'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\Block-1'
mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
load(fullfile(mafolder, mafilestruct.name)); % load SingleTargetKluverMAData
fs_ma = SingleTargetKluverMAData.SR;

% time indices for target onset, reach onset, touch screen, return and mouth
TargetTime = SingleTargetKluverMAData.TargetTime;
ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
MouthTimeix = SingleTargetKluverMAData.MouthTimeix;

timeixtbl_ma = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];

% extract the indices of the good trials
tag_goodreach = SingleTargetKluverMAData.goodix_reach;
tag_goodreturn = SingleTargetKluverMAData.goodix_return;
idx_goodtrials = find(tag_goodreach .* tag_goodreturn == 1);
timeixtbl_ma = timeixtbl_ma(idx_goodtrials,:);

tevent_ma = timeixtbl_ma{:,:}/fs_ma;
tevent_varNames = {'target', 'reach', 'touch', 'return', 'mouth'};

clear TargetTime ReachTimeix TouchTimeix ReturnTimeix MouthTimeix
end


function [tdur_ma, tdur_varNames] = readMAdur(onedaypath, block)
% read the MA duration time data from one file
%
%
%  Inputs:
%   onedaypath = 'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_013017'
%   block = 10  
%
%  Outputs:
%   tevent_ma : matrix of behavior time slots ntrials * 5
%   tevent_varNames: names of each col in tevent_ma 1 *5 cell 
%
mafolder = fullfile(onedaypath, ['Block-' num2str(block)]); %  'Y:\Animals2\Pinky\Recording\Processed\DataDatabase\Pinky_071417\Block-1'
mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
load(fullfile(mafolder, mafilestruct.name)); % load SingleTargetKluverMAData
tdur_ma = cat(2, SingleTargetKluverMAData.reaction_time, SingleTargetKluverMAData.reach_time, ...
    SingleTargetKluverMAData.manipulation_time,SingleTargetKluverMAData.return_time);


% extract the indices of the good trials
tag_goodreach = SingleTargetKluverMAData.goodix_reach;
tag_goodreturn = SingleTargetKluverMAData.goodix_return;
idx_goodtrials = find(tag_goodreach .* tag_goodreturn == 1);
%tdur_ma = tdur_ma(idx_goodtrials,:);

tdur_varNames = {'reaction_time', 'reach_time', 'manipulation_time', 'return_time'};

clear TargetTime ReachTimeix TouchTimeix ReturnTimeix MouthTimeix
end

