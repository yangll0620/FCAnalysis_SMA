function [ciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(ciCoh, T_chnsarea)
%
%   Input:
%       ciCoh: nchns * nchns * nf
%       T_chnsarea: nchns * 5
%
%
[nchns, ~, nf] = size(ciCoh);
chnPairNames = {};
ciCoh_flatten = zeros(nchns * (nchns -1)/2, nf);
ci = 0;
for chni = 1 : nchns -1
    for chnj = chni + 1  : nchns
        chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
        
        ci = ci + 1;
        ciCoh_flatten(ci, :) = ciCoh(chni, chnj, :);
    end
end