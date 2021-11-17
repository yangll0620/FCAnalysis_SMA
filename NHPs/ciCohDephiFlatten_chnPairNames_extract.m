function [ciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(ciCoh, deltaphis_allChnsTrials, T_chnsarea)
%
%   Input:
%       ciCoh: nchns * nchns * nf
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials(nsegs)
%       T_chnsarea: nchns * 5
%
%
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