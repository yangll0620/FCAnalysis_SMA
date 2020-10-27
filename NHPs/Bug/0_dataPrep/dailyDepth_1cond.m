function [t_areaDailyDepth]  = dailyDepth_1cond(xlsfile, sheetname, nCols, row_chn, row_area, row_iniDepth, row_recStart)
%
%
%   Inputs
%       xlsfile: the xls file (e.g /home/Bug/BugMasterDatabase.xlsx)
%       sheetname: the sheet name (e.g Depth of GM array_Moderate Chan)
%       nCols: the total column number, a scalar
%       row_chn:  the row number for chn 
%       row_area: the row number for area 
%       row_iniDepth: the row number for the initial depth
%       row_recStart: the row number for the recording start
%
%   Return:
%       t_areaDailyDepth: a table containing both area and daily depth (can be written directy using writetable)
%                         the variableNames are the date and channel numbers (e.g Date, chan4, chan19 et al)




% extract the chan* as variable names for all the following tables
varTypes = cell(1, nCols);
varTypes(:) = {'string'};
opts = spreadsheetImportOptions('NumVariables', nCols, ...
                                'VariableTypes', varTypes, ...
                                'DataRange', ['A' num2str(row_chn) ':' colNum2ExcelColName(nCols) num2str(row_chn)]);
varNames = readcell(xlsfile, opts, 'Sheet', sheetname);
for i = 1: length(varNames) % change from 9 to 'chan9'
    if isa(varNames{i},'double')
        varNames{i} = ['chan' num2str(varNames{i})];
    end
end
varNames = cellfun(@(x) strrep(x, ' ', ''), varNames, 'UniformOutput', false);
varNames{1} = 'Date';


% t_barea: the area for each channel
areaNames = readcell(xlsfile, 'FileType', 'spreadsheet', 'Sheet', sheetname, 'Range', ['A' num2str(row_area) ':' colNum2ExcelColName(nCols) num2str(row_area)]);
for i = 2: length(areaNames)
    if ismissing(areaNames{i})
        areaNames{i} = areaNames{i -1};
    end
end
areaNames{1} = '';
t_barea = cell2table(areaNames, 'VariableNames', varNames);
disp('channel number and brain area for each channel extracted.....')


% t_iniDepth: the initial depth for each channel
iniDepth = readmatrix(xlsfile, 'FileType', 'spreadsheet', 'Sheet', sheetname, 'Range', ['A' num2str(row_iniDepth) ':' colNum2ExcelColName(nCols) num2str(row_iniDepth)]);
t_iniDepth = array2table(iniDepth, 'VariableNames', varNames);




% t_dailyDepth: the daily depth table including Date as the first column
c_chnTypes = cell(1, nCols -1);
c_chnTypes(:) = {'double'};
varTypes = ['string', c_chnTypes];
opts = spreadsheetImportOptions('Sheet', sheetname, ... 
                                'VariableNames', varNames, ...
                                'VariableTypes', varTypes, ...
                                'DataRange', ['A' num2str(row_recStart)]);
t_records = readtable(xlsfile, opts);
rowmasks_depth = cellfun(@(x) ~isempty(x), t_records.Date);
t_recDepth = t_records(rowmasks_depth, :);
clear t_records rowmasks_depth c_chnTypes opts
disp('daily record extracted.....')



% fill depth for each channel
t_dailyDepth = t_recDepth(:, 'Date');
for i = 2 : length(varNames)
    varName = varNames{i};
    depth_1chn = t_recDepth{:, varName};
    iniDepth = t_iniDepth{:, varName};
    
    filled_depth = filled_dailyDepth(depth_1chn, iniDepth);
    
    t_dailyDepth = [t_dailyDepth table(filled_depth, 'VariableNames', {varName})];
    
    clear depth_1chn iniDepth filled_depth varName
end
disp('daily depth  filled .....')




% combine t_barea and t_dailyDepth

dailyDepthStr= string(t_dailyDepth{:, 2:end});
t_dailyDepthStr = array2table(dailyDepthStr, 'VariableNames', varNames(2:end));
t_dailyDepthStr = [t_dailyDepth(:, 1) t_dailyDepthStr];

t_areaDailyDepth = [t_barea; t_dailyDepthStr];
end




function depth_filled = filled_dailyDepth(depth_1chn, iniDepth)
% filled the depth, if empty, use the previous one, 
%                   otherwise depth_1chn(i) + initDepth
% 
% Input
%       depth_1chn: the column Number ndays * 1
%       iniDepth: initial depth, scale
% output:
%       depth_filled: filled depth vector ndays * 1

depth_filled = zeros(size(depth_1chn));

% fill the first day
i = 1;
if isnan(depth_1chn(i))
    depth_filled(i) = iniDepth;
else
    depth_filled(i) = iniDepth + depth_1chn(i);
end

for i = 2: length(depth_1chn)
    
    if isnan(depth_1chn(i))
        depth_filled(i) = depth_filled(i-1);
        
    else
        depth_filled(i) = iniDepth + depth_1chn(i);
    end
end
end


function cCol= colNum2ExcelColName(colNum)
% colnum to excel column name
% Exp: 1-> 'A',  27 -> 'AA'
% Input
%       colNum: the column Number
% output:
%       cCol: the transformed excel column name

    div = floor(colNum/26);
    re = mod(colNum, 26);


    cRe = char(re + 'A' - 1);

    if div ==  0
        cDiv = '';
    else
       cDiv = char(div + 'A' - 1) ;
    end
    cCol = [cDiv cRe];
end












