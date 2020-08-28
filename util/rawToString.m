function [ cellArray, stringArray ] = rawToString( cellArrayIn )
%RAW_TO_STRING Summary of this function goes here
%   Detailed explanation goes here
 % changing format to agree with raw data from function readxls
    cellArray = cellArrayIn;
    nrows = size(cellArray,1);
    colMax = 0;
    for row=1:nrows
        ncols = size(cellArray,2);
        if ncols > colMax
            colMax=ncols;
        end
    end
    
    for row=1:nrows
        ncols = size(cellArray,2);
        if ncols < colMax
            cellArray(row,colMax)={[]};
        end
    end
 
    for row=1:nrows
        for col=1:ncols
            thisCell =cellArray(row,col);
            if isempty(thisCell{1})
              
                cellArray(row,col)={[NaN]};
                stringArray(row,col)={''};
            else
                stringArray(row,col) = thisCell;
                if ~isempty( str2num(thisCell{1}) )
                    cellArray(row,col)={str2num(thisCell{1})}; 
                end
            end
        end
    end

end

