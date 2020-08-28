function m5_restData_noM1PMC_avgEachArea()
%   extract avg lfp data in each area (DBS remained all contacts)
%
%   Processing steps as follows:
%       1. for M1 and PMC, extract the useful channels for each exp date
%
%       2. average lfp for each area



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

% input folder: extracted raw STK data
inputfolder = fullfile(codecorresParentfolder, 'm4_restData_segment_chnArea');


% depth and reference channel for M1
depth_M1Layer5 = [1.25 1.75] * 8;
recordchnref_M1 = 58;

% depth and reference channel for PMC
depth_PMCLayer5 = [1.25 1.75] * 8;
recordchnref_PMC = 50;


%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'noM1PMC_avgArea';




%% Code Start Here
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
close all;
f = waitbar(0, ['Processing lfp data...']);

for i = 1 : nfiles
    % wait bar
    waitbar(i/nfiles,f,['Processing  file ' num2str(i) '/' num2str(nfiles)]);
    
    % load lfpdata: nchns * ntemp * ntrials
    filename = files(i).name;
    load(fullfile(inputfolder, filename), 'lfpdata', 'fs', 'T_chnsarea');
    
    
    % extract the date of the exp
    idx = strfind(filename, '_bktdt');
    dateofexp = datetime(filename(idx-8: idx-1), 'InputFormat', 'yyyyMMdd');
    
    
    %%% ----- extract averaged lfp within each brain area except M1, PMC and DBS ------%%%
    
    % extract the unique GM Areas, excluding the empty, combination (e.g.
    % 'VA/Vlo/STN' ), 'M1' , 'PMC'  and the DBS
    
    GMChnAreas = T_chnsarea.brainarea;
    uniqGMAreas = unique(T_chnsarea.brainarea);
    idxs = cellfun(@(x) isempty(x)||contains(x,'M1')||contains(x,'PMC')||contains(x,'/')||contains(x,'STN')||contains(x,'GP'), uniqGMAreas);
    uniqGMAreas(idxs) = [];
    
    
    % average lfp across each remaining GM area (except M1 and PMC) in uniqGMAreas
    newrows_GM = T_chnsarea(1,:);
    
    avglfp_remainGM = [];
    for i = 1: length(uniqGMAreas)
        GMArea = uniqGMAreas(i);
        
        % find the channel numbers of GMArea in GMChnAreas(i.e. in lfpsegs_GM)
        indices = find(cellfun(@(x) strcmp(x, GMArea), GMChnAreas));
        
        % average across GMArea
        avglfp = mean(lfpdata(:,indices, :),2);
        
        % concatenate to avglfpsegs_GM
        avglfp_remainGM = cat(2, avglfp_remainGM, avglfp);
        
        if i == 1
            newrows_GM.chni(1) = i;
            newrows_GM.brainarea{1} = GMArea;
            newrows_GM.recordingchn(1) = nan;
            newrows_GM.electype{1} = 'Gray Matter';
            newrows_GM.notes{1} = indices;
        else
            newrows_GM.chni(end+1) = i + 2;
            newrows_GM.brainarea{end} = GMArea;
            newrows_GM.recordingchn(end) = nan;
            newrows_GM.electype{end} = 'Gray Matter';
            newrows_GM.notes{end} = indices;
            
        end
        
        
        clear GMArea indices avglfpsegs
    end
    avglfp_GM = cat(2, avglfp_remainGM);
    
    
    
    %%% ----- deal with lfpdata and T_chnsarea  ------%%%
    
    rowidx_DBS = (T_chnsarea.brainarea == "STN" | T_chnsarea.brainarea == "GP");
    
    lfp_DBS = lfpdata(:, rowidx_DBS, :);
    lfpdata = cat(2, avglfp_GM, lfp_DBS);
    
    
    rows_DBS = T_chnsarea(rowidx_DBS,:);
    T_chnsarea = [newrows_GM;rows_DBS];
    
    
    % save
    savefilename = [animal '_' filename(4:end)]; 
    
    save(fullfile(savefolder, savefilename), 'lfpdata','fs', 'T_chnsarea');
    
    clear newrows_M1PMC
    clear lfpM1 lfpPMC recordchnUsed_PMC recordchnUsed_M1 rowidx_M1PMC
    clear lfpdata fs T_chnsarea T_idxevent
end



function recchns_useful = useful_recordchn_extract(brainarea, depth_useful, dateofexp)
% extract the recording channel number in brain area in the range depth_useful in the date dateofexp
%
%     args:
%           brainarea: the brain area, i.e. 'M1', 'SMA'
%
%           depth_useful (vector): the depth range of the useful channel, i.e. depth_useful = [1.25 1.75] * 8
%
%           dateofexp (datetime): the exp date, i.e. dateofexp = datetime('04/02/2019', 'Format','MM/dd/yyyy')
%
%    return:
%           recchns_useful: the useful recording channels, a vector. (e.g 93 or [90 91])



%% start here

% find the codefolder
filepath_currCode = mfilename('fullpath');
idx = strfind(filepath_currCode, 'code');
codefolder = filepath_currCode(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

cond = parsePDCondition(datenum(dateofexp), 'Bug');

% channel depth information for each exp date
[datafolder, ~, ~, ~] = exp_subfolders();
if strcmp(cond, 'mild')
    filename_chanDepth = 'Bug_channelDepth_mild.csv';
else
    if strcmp(cond, 'normal')
        filename_chanDepth = 'Bug_channelDepth_normal.csv';
    end
end
file_chanDepth = fullfile(datafolder, 'Bug', filename_chanDepth);

% return the import options
opts = detectImportOptions(file_chanDepth);

varnames = opts.VariableNames;
nvars = length(varnames);

% read format for each column
fmt = ['%{MM/dd/yy}D' repmat('%f', 1, nvars-1)];

% read the channel depth information
warning off
T_chanDepth = readtable(file_chanDepth, 'Format',fmt);
warning on

% extract the column indices for M1
idx_col = cellfun(@(x) contains(x, brainarea), varnames);


% extract the row index for the exp date
idx_row = find(T_chanDepth{:, 1} == dateofexp);
if isempty(idx_row)
    disp([datestr(dateofexp) ': idx_row is empty']);
    
    recchns_useful = [];
    return;
end

if length(idx_row)>1
    disp([datestr(dateofexp) ': idx_row has more than 1 row!' ]);
    
    recchns_useful = [];
    return;
end


% extract the useful M1 channels e.g. [1.25 1.75]mm for the exp date
T_area_dateofexp = [T_chanDepth(1, idx_col); T_chanDepth(idx_row, idx_col)];
idx_useful = find(T_area_dateofexp{2, :}>=depth_useful(1) & T_area_dateofexp{2, :} <=depth_useful(2));
T_useful = T_area_dateofexp(:,idx_useful);

% find the correponding lfp data for the particular M1 channel in the date of exp
recchns_useful = [];
for i = 1: width(T_useful)
    recchns_useful = [recchns_useful; T_useful{1, i}];
end








