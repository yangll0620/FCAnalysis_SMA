function [convertedMatrix, varNamesCell]= conv_tbl2matrix(T)
%% convert table object into cell and varNamesCell as python could not read the matlab table objects
% 
%   Input:
%       T: matlab table object, the value of T should be a matrix
%
%	Output:
%		convertedMatrix: converted matrix, each column contains each column value of T
%		varNamesCell:  cells in which each cell contains each column variablename

convertedMatrix = T{:,:};
varNamesCell = T.Properties.VariableNames;