%% Usage example
%  modify the input parameters accordingly
%       xls master file 
%       save related information
%       Setup for normal, mild and moderate
%
%
%% input parameters
file_xlsmaster = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/data/Bug/BugMasterDatabase.xlsx';

% save
savefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/NHPs/Bug/0_dataPrep';
savename_prefix = 'Bug_dailyDepth';
savefile = fullfile(savefolder, [savename_prefix '.xlsx']);

% setup for normal
sheet_normal = 'Depth of GM array_Normal Channe';
nCols_normal = 36; % the total column number
row_area_normal = 19; row_chn_normal = 20;   % the row number for recording area, and channel number
row_iniDepth_normal= 21; row_recStart_normal = 22; % the row number for recording initial depth, and recording start


% setup for mild
sheet_mild = 'Depth of GM array_Mild Channels';
nCols_mild = 93; % the total column number
row_area_mild = 2;     row_chn_mild = 3;   % the row number for recording area, and channel number
row_iniDepth_mild= 5;  row_recStart_mild = 8; % the row number for recording initial depth, and recording start


% setup for moderate
sheet_moderate = 'Depth of GM array_Moderate Chan';
nCols_moderate = 94; % the total column number
row_area_moderate = 1;          row_chn_moderate = 2;   % the row number for recording area, and channel number
row_iniDepth_moderate= 8;   row_recStart_moderate = 10; % the row number for recording initial depth, and recording start




%% coding  running here
% normal case
disp(' ..... Dealing Normal .........')
t_areaDailyDepth_normal = dailyDepth_1cond(file_xlsmaster, sheet_normal, nCols_normal, row_chn_normal, row_area_normal, row_iniDepth_normal, row_recStart_normal);
writetable(t_areaDailyDepth_normal, savefile, 'sheet', 'normal');





% mild case
disp(' ..... Dealing mild .........')
t_areaDailyDepth_mild = dailyDepth_1cond(file_xlsmaster, sheet_mild, nCols_mild, row_chn_mild, row_area_mild, row_iniDepth_mild, row_recStart_mild);
writetable(t_areaDailyDepth_mild, savefile, 'sheet', 'mild');


% moderate case
disp(' ..... Dealing moderate .........')
t_areaDailyDepth_moderate = dailyDepth_1cond(file_xlsmaster, sheet_moderate, nCols_moderate, row_chn_moderate, row_area_moderate, row_iniDepth_moderate, row_recStart_moderate);
writetable(t_areaDailyDepth_moderate, savefile, 'sheet', 'moderate');


disp(['daily depth stored in ' savefile ', with sheetname normal, mild and moderate'])
