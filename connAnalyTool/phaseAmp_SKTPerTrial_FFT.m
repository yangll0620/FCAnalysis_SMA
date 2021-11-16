function [phis, amps, f_selected] = phaseAmp_SKTPerTrial_FFT(lfptrials, fs, f_AOI)
%
%   ref: https://yll0620.medium.com/functional-connectivity-measurement-16423fee3581
%      return the phase and amp for each trial data
%
%   Input
%       lfpdata:  ntemp * ntrials
%
%       f_AOI: frequencies duration of interest, i.e. f_AOI = [8 40]
%
%   Outputs:
%       phis: phase angle in radians for each trial (nf * ntrials)
%
%       amps: magnitude for each trial (nf * ntrials)
%
%       f_selected: selected frequencies (nf * 1)


[ntemp, ntrials] = size(lfptrials);
phis = [];
amps = [];
for triali = 1: ntrials
    x = lfptrials(: , triali);
    
    [Sx, fx, ~, ~] = spectrogram(x, ntemp, 0,[],fs); % Sx: nf * 1
    
    if triali == 1
        freqs = fx;
        idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
        f_selected = round(freqs(idx_f), 2);
        
        clear freqs
    end
    
    phi = angle(Sx(idx_f, :)); % phix : nf *1
    amp = abs(Sx(idx_f, :)); % amp: nf * 1
    
    phis = cat(2, phis, phi);
    amps = cat(2, amps, amp);
    
    clear x Sx fx phi amp
end
