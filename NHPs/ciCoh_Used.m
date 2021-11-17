function [ciCoh_flatten_used, chnPairNames_used]= ciCoh_Used(chnPairNames, ciCoh_flatten, removed_chns)
%
% Input:
%   ciCoh_flatten: npairs * nf
%   

removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
clear M1DBS_mask STN2GP_mask

ciCoh_flatten_used = ciCoh_flatten(usedChnPairsMask, :);
chnPairNames_used = chnPairNames(usedChnPairsMask);