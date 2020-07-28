function [psd, ts, F] = spectrogram_calc(lfp_phase, twin, toverlap, fs, t0)
%% spectrogram 
% 
%   inputs:
%       lfp_phase:  ntemp * ntrials
%
%       t0: the corresponding time slot for the first temporal sample
%
%   return:
%        psd: nfs * nts * ntrials

nwin = round(twin * fs);
noverlap = round(toverlap * fs);
psd = [];
for triali = 1 : size(lfp_phase, 2)
    lfp_1trial = lfp_phase(:,triali);
    
    [~, F, T, P] = spectrogram(lfp_1trial,nwin,noverlap,nwin,fs);
    
    
    % psd: nfs * nts * ntrials
    psd = cat(3, psd, P);
    
    clear lfp_1trial P
end

% time in ms
ts = T - t0;