function [iCoh, f_selected] = imCohRest_FFT_NormalizedAMP(lfpsegdata, twin, toverlap, fs, f_AOI)
% 
%   Input
%       lfpdata: chns * ntemp * nsegs
%
%       f_AOI: frequencies of interest, i.e. f_AOI = [8 40]
%       
%       t_AOI: time duration used to extract mean iCOH
%
%       fs: sample rate
%
%       twin, toverlap: time duration to calculate ciCOH, i.e. twin = 0.2; toverlap = 0.15;

nsegs = size(lfpsegdata, 3);
for segi = 1 : nsegs
    lfpdata = lfpsegdata(:, :, segi);
    nchns = size(lfpdata, 1);
    
    for chni = 1 : nchns-1
        lfpi = squeeze(lfpdata(chni, :));
        [Si, fi, ~, ~] = spectrogram(lfpi, round(twin * fs), round(toverlap * fs),[],fs); % Si: nf * nt
        
        % extract idx_f, f_selected and iCoh_eachSeg at the first time
        if segi == 1 && chni == 1
            idx_f = (fi>=f_AOI(1) &  fi<=f_AOI(2));
            f_selected =  round(fi(idx_f), 2);
            
            nf = length(f_selected);
            crossDensity_eachSeg = zeros(nchns, nchns, nf, nsegs);
            
            clear nf
        end
        
        phii = angle(Si(idx_f, :));
        ampi = abs(Si(idx_f, :));
        normAMPi = mean((ampi .* ampi), 2);
        
        for chnj = chni + 1 : nchns
            lfpj = squeeze(lfpdata(chnj, :));
            
            [Sj, ~, ~, ~] = spectrogram(lfpj, round(twin * fs), round(toverlap * fs),[],fs); % Sj: nf * nt
            phij = angle(Sj(idx_f, :));
            ampj = abs(Sj(idx_f, :));
            normAMPj = mean((ampj .* ampj), 2);
            
            crossDensity_eachSeg(chni, chnj, :, segi) = mean(ampi .* ampj .* exp(1i * (phii - phij)), 2) ./ sqrt(normAMPi .* normAMPj);
            crossDensity_eachSeg(chnj, chni, :, segi) = crossDensity_eachSeg(chni, chnj, :, segi);
            clear lfpj Sj phij ampj
        end
        
        clear lfpi Si fi  phii ampi normAMPi
    end
    
    clear lfpdata nchns
end
iCoh = abs(imag(mean(crossDensity_eachSeg, 4))); % iCoh: nchns * nchns * nf