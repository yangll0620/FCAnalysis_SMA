function m1_dailyDepth()
%% Extract the daily depth from xlsmaster sheet for normal, mild and moderate
%  modify the input parameters accordingly
%       xls master file 
%       save related information
%       Setup for normal, mild and moderate
%
%

%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');


% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

% datafolder
[datafolder, ~, ~, ~] = exp_subfolders();


% animal
[fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
animal = codecorresfolder(fi + length('NHPs/'):j);


%% input parameters
file_xlsmaster = fullfile(datafolder, animal, 'BugMasterDatabase.xlsx');


% save
savefolder = codecorresfolder;
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



% remove the sheet 1, 2 and 3
book = matlab.io.spreadsheet.internal.createWorkbook('xlsx', savefile, false);
book.removeSheet('Sheet1')
book.removeSheet('Sheet2')
book.removeSheet('Sheet3')
book.save(savefile)
disp(['daily depth stored in ' savefile ', with sheetname normal, mild and moderate'])
