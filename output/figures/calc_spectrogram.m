function [psd_allchns, freqs, times] = calc_spectrogram(lfp_phase_trials, fs, tdur_trial)
%
% Inputs:
%    lfp_phase_trials: nchns * ntemp * ntrials
%
% Return:
%   psd_allchns: nf * nt * nchns
%   freqs: nf * 1
%   times: 1 * nt
  
twin = 0.2;
toverlap = 0.18;
f_AOI = [8 40];
t_AOI = [-0.5 0.5];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);

psd_allchns = [];
for chi = 1 : size(lfp_phase_trials, 1)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfp_phase_trials, 3)
        x = lfp_phase_trials(chi, :, tri);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end
end