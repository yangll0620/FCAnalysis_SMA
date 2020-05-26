function [ lfp ] = get_lfp_removeBadChns( lfp_array,  m1array_to_remove)
% get_lfp_removeBadChns Get LFP from a 10 by 10 array of lfp_array
% modified from Ying's code: get_lfp_mean(lfp_array, m1array_to_remove)  
%
%   Args:
%       lfp_array: 1 * 10 cell, lfp_array{rowi}{colj} (ntemp * 1) is time series lfp data of
%                  one channel 
%
%       m1Array_to_remove: nchns_removed * 2 (1: rowi,  2: colj)
% 
%   Return:
%       lfp: ntemp * nchns, nchns = 100 - length(m1_array_to_remove)

   
lfp = [];
for row=1:10
    for col=1:10
        if isempty (strmatch([row col], m1array_to_remove))
            lfp = cat(2,lfp, lfp_array{row}{col});
        end
    end
end
