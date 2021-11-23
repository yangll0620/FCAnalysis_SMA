function [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, ciCoh_flatten, deltaphis_flatten, removed_chns, varargin)
%
% Input:
%   ciCoh_flatten: npairs * nf
%   deltaphis_flatten: npairs * nf * ntrials
%
%       Name-Value: 
%           'codesavefolder' - code saved folder


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
clear M1DBS_mask STN2GP_mask

ciCoh_flatten_used = ciCoh_flatten(usedChnPairsMask, :);
deltaphis_flatten_used = deltaphis_flatten(usedChnPairsMask, :, :);
chnPairNames_used = chnPairNames(usedChnPairsMask);