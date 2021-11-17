function [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh)
%
%   1. fit a normal distribution based on psedociCohs
%
%   2. extract pvalues using permutation test and extract significate ciCoh using Benjamini & Hochberg (1995) procedure 
%
%   Input:
%       psedociCohs: nchns * nchns * nf * nshuffles
%
%   Output:
%       sigciCoh: nchns * nchns * nf (if not sig, set 0; otherwise remain)

% fit a normal distribution
[nchns, ~, nf, ~] = size(psedociCohs);
mus = zeros(nchns, nchns, nf);
stds = zeros(nchns, nchns, nf);
for fi = 1: nf
    for chni = 1 : nchns -1
        for chnj = chni + 1 : nchns
            x = squeeze(psedociCohs(chni, chnj, nf, :));
            pd = fitdist(x,'Normal');
            mus(chni, chnj, fi) = pd.mu;
            stds(chni, chnj, fi) = pd.std;
            clear x pd
        end
    end
end

% pvalues using permutation test
[nchns, ~, nf] = size(ciCoh);
pvals = zeros(size(ciCoh));
for fi = 1 : nf
    for chni = 1: nchns -1
        for chnj = chni : nchns
            mu1 = mus(chni, chnj, fi);
            std1 = stds(chni, chnj, fi);
            pd = makedist('Normal','mu',mu1,'sigma',std1);
            
            x = ciCoh(chni, chnj, fi);
            pvals(chni, chnj, fi) = (1 - cdf(pd,abs(x)));
            
            clear x
            clear mu1 std1 pd
        end
    end
end
% Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR)
[h, ~, ~, ~]=fdr_bh(pvals);

% set values not significant as 0
sigciCoh = ciCoh;
sigciCoh(h==0) = 0;