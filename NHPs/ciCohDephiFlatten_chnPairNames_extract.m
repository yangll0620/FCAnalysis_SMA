function [ciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(ciCoh, deltaphis_allChnsTrials, T_chnsarea, varargin)
%
%   Input:
%       ciCoh: nchns * nchns * nf
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials(nsegs)
%       T_chnsarea: nchns * 5
%
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


[nchns, ~, nf] = size(ciCoh);
chnPairNames = {};
ciCoh_flatten = zeros(nchns * (nchns -1)/2, nf);
ntrials = size(deltaphis_allChnsTrials, 4);
deltaphis_flatten = zeros(nchns * (nchns -1)/2, nf, ntrials);
ci = 0;
for chni = 1 : nchns -1
    for chnj = chni + 1  : nchns
        chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
        
        ci = ci + 1;
        ciCoh_flatten(ci, :) = ciCoh(chni, chnj, :);
        deltaphis_flatten(ci, :, :) = deltaphis_allChnsTrials(chni, chnj, :, :);
    end
end