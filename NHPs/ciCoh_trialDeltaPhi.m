function [deltaphis_allChnsTrials, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfptrials, fs, f_AOI)
%
%
%   
%   Output:
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials
%       ciCoh:  nchns * nchns * nf
%       f_selected: nf * 1



% extract ciCoh
[ciCoh, f_selected] = ciCohSKTAllchns_FFT_NoAmp(lfptrials, fs, f_AOI);

% extract deltaphis for all chns and trials
[nchns, ~, ntrials] = size(lfptrials);
nf = length(f_selected);
deltaphis_allChnsTrials = zeros(nchns, nchns, nf, ntrials);
for chni = 1: nchns - 1
    lfptriali = squeeze(lfptrials(chni, :, :));
    [phisi, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfptriali, fs, f_AOI);
    for chnj = chni + 1 : nchns
        lfptrialj = squeeze(lfptrials(chnj, :, :));
        [phisj, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfptrialj, fs, f_AOI);
        deltaphis_allChnsTrials(chni, chnj, :, :) = phisi - phisj;
        clear lfptrialj phisj
    end
    clear lfptriali phisi
end