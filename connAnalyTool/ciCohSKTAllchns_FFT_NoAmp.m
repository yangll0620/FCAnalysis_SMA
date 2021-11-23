function [ciCoh_NoAmp, f_selected] = ciCohSKTAllchns_FFT_NoAmp(lfpdata, fs, f_AOI)
%
%   ref: https://yll0620.medium.com/functional-connectivity-measurement-16423fee3581 
%      treat ntemp as stationary time series, i.e. twin = ntemp, toverlap = 0
%
%   Input
%       lfpdata: nchns * ntemp * ntrials
%
%       f_AOI: frequencies duration of interest, i.e. f_AOI = [8 40]
%
%   Outputs:
%       ciCoh_NoAmp: absolute imaginery coherence nchns * nchns * nf, Upper triangle


nchns = size(lfpdata, 1);

% extract phis and amps for each trial and each channel
for chni = 1 : nchns
    lfptrials = squeeze(lfpdata(chni, :, :));
    [phis, ~, f_selected] = phaseAmp_SKTPerTrial_FFT(lfptrials, fs, f_AOI);
    
    if chni == 1
        [nf, ntrials] = size(phis);
        phis_allchns = zeros(nchns, nf, ntrials); % phis_allchns: nchns * nf * ntrials
        clear nf ntrials
    end
    
    phis_allchns(chni, :, :) = phis;
    
    clear lfptrials phis  
end


% calculate ciCoh
nf = size(phis_allchns, 2);
ciCoh_NoAmp = zeros(nchns, nchns, nf);
for chni = 1: nchns -1
    phii = squeeze(phis_allchns(chni, :, :));
    for chnj = chni + 1 : nchns
        phij = squeeze(phis_allchns(chnj, :, :));
        
        deltaphi = phii - phij;
        coh = mean(exp(1i * deltaphi), 2);
        ciCoh_NoAmp(chni, chnj, :) = abs(imag(coh) ./ sqrt((1 - real(coh).^2)));

        clear phij deltaphi
    end
    clear phii ampi Gii
end