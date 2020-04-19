function [convertedCell, varNamesCell]= conv_tbl2Cell(T)
%% convert table object into cell and varNamesCell as python could not read the matlab table objects
% 
%   Input:
%       T: matlab table object
%
%	Output:
%		convertedCell: converted cells, each cell member contains each column value
%		varNamesCell:  cells in which each cell contains each column variablename

convertedCell = {};
varNamesCell = {};
for widi = 1 : width(T)
    
    varNamesCell = [varNamesCell;T.Properties.VariableNames{widi}];
    convertedCell = [convertedCell;{T{:,widi}}];
end
