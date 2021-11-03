function [iCoh_acrossTime, f_selected] = imCohSKT_FFT_NormalizedAMP(lfpdata, twin, toverlap, fs, f_AOI, t_AOI, tdur_trial)
% 
%   Input
%       lfpdata: chns * ntemp * nsegs
%
%       f_AOI, t_AOI: frequencies and time duration of interest, i.e. f_AOI
%       = [8 40], t_AOI = [-0.2 0]
%
%       twin, toverlap: time duration to calculate ciCOH, i.e. twin = 0.2; toverlap = 0.15;
%
%       tdur_trial: time duration for used trial data, i.e tdur_trial = [-0.5 0.5]
%
%   Outputs:
%       iCoh_acrossTime: imaginery coherence nchns * nchns * nf



[nchns, ~, ntrials] = size(lfpdata);
for chni = 1 : nchns-1
    lfptrialsi = squeeze(lfpdata(chni, :, :));
    for chnj = chni : nchns
        lfptrialsj = squeeze(lfpdata(chnj, :, :));
        
        cross_density_sum = 0;
        densityX_sum = 0;
        densityY_sum = 0;
        for triali = 1: ntrials
            x = lfptrialsi(: , triali);
            y = lfptrialsj(: , triali);
            
            [Sx, fx, tx, ~] = spectrogram(x, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
            [Sy, ~, ~, ~] = spectrogram(y, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
            
            if chni == 1 && chnj == chni && triali == 1
                freqs = fx;
                times = tx + tdur_trial(1);
                idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
                idx_t = (times>=t_AOI(1) &  times<=t_AOI(2));
                f_selected = round(freqs(idx_f), 2);
                clear freqs times
            end
            
            phix = angle(Sx(idx_f, idx_t));
            phiy = angle(Sy(idx_f, idx_t));
            ampx = abs(Sx(idx_f, idx_t));
            ampy = abs(Sy(idx_f, idx_t));
            cross_density_sum = cross_density_sum + ampx .* ampy .* exp(1i * (phix - phiy));
            densityX_sum = densityX_sum + ampx .* ampx;
            densityY_sum = densityY_sum + ampy .* ampy;
            
            clear x y Sx fx tx Sy
            clear phix phiy
            clear ampx ampy
        end
        
        if chni == 1 && chnj == chni + 1
            [nf, nt] = size(cross_density_sum);
            cross_densitys = zeros(nchns, nchns, nf, nt);
            densitysX = zeros(nchns, nchns, nf, nt);
            densitysY = zeros(nchns, nchns, nf, nt);
            clear nf nt
        end
        cross_densitys(chni, chnj, :, :) = cross_density_sum / ntrials; % cross_densitys: nchns * nchns * nf  * nt
        cross_densitys(chnj, chni, :, :) = cross_densitys(chni, chnj, :, :);
        
        densitysX(chni, chnj, :, :) = densityX_sum / ntrials; % cross_densitys: nchns * nchns * nf  * nt
        densitysX(chnj, chni, :, :) = densitysX(chni, chnj, :, :);
        
        densitysY(chni, chnj, :, :) = densityY_sum / ntrials; % cross_densitys: nchns * nchns * nf  * nt
        densitysY(chnj, chni, :, :) = densitysY(chni, chnj, :, :);
        
        
        clear cross_density_sum lfptrialsj
        clear densityX_sum densityY_sum
    end
    
    clear lfptrialsi
end
iCoh = imag(cross_densitys ./ sqrt(densitysX .* densitysY));
iCoh_acrossTime = abs(mean(iCoh, 4));