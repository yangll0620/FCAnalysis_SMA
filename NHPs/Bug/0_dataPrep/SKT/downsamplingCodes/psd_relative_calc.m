function psd_rel = psd_relative_calc(psd, psd_base)
% relative psd: 10log10(psd/psd_base)
%
%   parameters:
%       psd: nfs * nts_phase ** ntrials
%       
%       psd_base: nfs * nts_base * ntrials
% 
%   return:
%       pad_rel: nts * nfs * ntrials


% base psd average across the time
psd_base = mean(psd_base, 2);

nts = size(psd, 2);

% relative psd
psd_rel = 10 * log10(psd ./ repmat(psd_base, 1, nts, 1));