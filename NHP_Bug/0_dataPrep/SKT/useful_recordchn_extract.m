function recordingchns_useful = useful_recordchn_extract(brainarea, depth_usefu, dateofexp)
% extract the recording channel number in brain area in the range depth_useful in the date dateofexp 
%
%     args:
%           brainarea: the brain area, i.e. 'M1', 'SMA'
%
%           depth_useful (vector): the depth range of the useful channel, i.e. depth_useful = [1.25 1.75] * 8
%
%           dateofexp (datetime): the exp date, i.e. dateofexp = datetime('04/02/2019', 'Format','MM/dd/yyyy')
%  



%% start here

addpath('/Users/linglingyang/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code/util');

cond = parsePDCondition(datenum(dateofexp), 'Bug');

% channel depth information for each exp date
datafolder = '/Users/linglingyang/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/data';
if strcmp(cond, 'mild')
    filename_chanDepth = 'Bug_channelDepth_mild.csv';
else
    if strcmp(cond, 'normal')
        filename_chanDepth = 'Bug_channelDepth_normal.csv';
    end
end    
file_chanDepth = fullfile(datafolder, filename_chanDepth);

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
    
    recordingchns_useful = [];
    return;
end

if length(idx_row)>1
    disp([datestr(dateofexp) ': idx_row has more than 1 row!' ]);
    
    recordingchns_useful = [];
    return;
end


% extract the useful M1 channels e.g. [1.25 1.75]mm for the exp date 
T_area_dateofexp = [T_chanDepth(1, idx_col); T_chanDepth(idx_row, idx_col)];
idx_useful = find(T_area_dateofexp{2, :}>=depth_usefu(1) & T_area_dateofexp{2, :} <=depth_usefu(2));
T_useful = T_area_dateofexp(:,idx_useful);

% find the correponding lfp data for the particular M1 channel in the date of exp
recordingchns_useful = [];
for i = 1: width(T_useful)
    recordingchns_useful = [recordingchns_useful; T_useful{1, i}];
end




