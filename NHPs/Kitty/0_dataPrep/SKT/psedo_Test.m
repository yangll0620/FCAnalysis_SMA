function [mus, stds] = psedo_Test(lfp1, lfp2, shuffleN, fs, twin, toverlap, f_AOI)
% lfp1, lfp2: ntemp * ntrials
%
% return:
%       mus, stds: the mu and std of the fitted normal distribution (nf * nt)

ntrials = size(lfp2, 2);
lfp2_origin = lfp2;

for si = 1 : shuffleN
    % shuffle lfp2
    idx = randperm(ntrials);
    lfp2 = lfp2_origin(:, idx); 
    
    cross_density_sum = 0;
    for triali = 1: ntrials
        x = lfp1(: , triali);
        y = lfp2(: , triali);
        
        [Sx, fx, tx, ~] = spectrogram(x, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
        [Sy, fy, ty, ~] = spectrogram(y, round(twin * fs), round(toverlap * fs),[],fs); % Sy: nf * nt
        
        if triali == 1
            freqs = fx;
            times = tx;
            idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
            f_selected =  freqs(idx_f);
        end
        
        phix = angle(Sx(idx_f, :));
        phiy = angle(Sy(idx_f, :));
        cross_density_sum = cross_density_sum + exp(1i * (phix - phiy));
        
        clear x y Sx fx tx Sy fy ty
        clear phix phiy
    end
    
    if si == 1
        [nf, nt] = size(cross_density_sum);
        cross_densitys = zeros(shuffleN, nf, nt);
        clear nf nt
    end
    
    cross_densitys(si, :, :) = cross_density_sum / ntrials;
    
    clear idx cross_density_sum
end

% abs imaginary of Coherency
psedo_iCoh = imag(cross_densitys);

% fit a normal distribution to psedo_iCoh
[~, nf, nt] = size(psedo_iCoh);
mus = zeros(nf,nt);
stds = zeros(nf, nt);
for fi = 1 : nf
    for ti = 1 : nt
        pd = fitdist(squeeze(psedo_iCoh(:, fi,ti)),'Normal');
        mus(fi, ti) = pd.mu;
        stds(fi, ti) = pd.sigma;
        
        clear pd
    end
end