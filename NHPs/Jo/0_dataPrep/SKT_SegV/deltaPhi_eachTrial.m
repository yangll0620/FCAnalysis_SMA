function [deltaPhis, f_selected] = deltaPhi_eachTrial(lfptrialsi, lfptrialsj, fs, twin, toverlap, f_AOI, t_AOI, tdur_trial)
%
%   Input:
%       lfptrialsi, lfptrialsj: lfp data for channel i and j (ntemp * ntrials)
%
%       f_AOI, t_AOI: frequencies and time duration of interest, i.e. f_AOI
%       = [8 40], t_AOI = [-0.2 0]
%
%       twin, toverlap: time duration to calculate ciCOH, i.e. twin = 0.2; toverlap = 0.15;
%
%       tdur_trial: time duration for used trial data, i.e tdur_trial = [-0.5 0.5]
%
%   Outputs:
%       deltaPhis: deltaPhi for each trial (nfs  * ntrials)

deltaPhis = [];
[ntemp, ntrials] = size(lfptrialsi);
nwin = ntemp; % nwin = twin * fs;
noverlap = (0 * fs);  % noverlap = (toverlap * fs);
for triali = 1: ntrials
    x = lfptrialsi(: , triali);
    y = lfptrialsj(: , triali);
    
    [Sx, fx, tx, ~] = spectrogram(x, nwin, noverlap,[],fs); % Sx: nf * nt
    [Sy, ~, ~, ~] = spectrogram(y, nwin, noverlap,[],fs); % Sy: nf * nt
    
    if triali == 1
        freqs = fx;
        times = tx + tdur_trial(1);
        idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
        idx_t = (times>=t_AOI(1) &  times<=t_AOI(2));
        f_selected = round(freqs(idx_f), 2);
        clear freqs times
    end
    
    phix = angle(Sx(idx_f, idx_t));
    phiy = angle(Sy(idx_f, idx_t));
    deltaPhi = mean(phix - phiy, 2);
    deltaPhis = cat(2, deltaPhis, deltaPhi);
    
    
    clear x y Sx fx tx Sy
    clear phix phiy deltaPhi
    clear ampx ampy
end