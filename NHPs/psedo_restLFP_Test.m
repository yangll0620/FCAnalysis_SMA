function [mus, stds] = psedo_restLFP_Test(lfp1, lfp2, shuffleN, fs, twin, toverlap, f_AOI)
% lfp1, lfp2: 1 * ntemp
%
% return:
%       mus, stds: the mu and std of the fitted normal distribution (nf * 1)

ntemp = size(lfp2, 2);
lfp2_origin = lfp2;

for si = 1 : shuffleN
    
    % shuffle lfp2 along ntemp
    idx = randperm(ntemp);
    lfp2 = lfp2_origin(:, idx);
    
    [S1, fi, ~, ~] = spectrogram(lfp1, round(twin * fs), round(toverlap * fs),[],fs); % Si: nf * nt
    [S2, ~, ~, ~] = spectrogram(lfp2, round(twin * fs), round(toverlap * fs),[],fs); % Sj: nf * nt
    
    % extract idx_f, f_selected and iCoh_eachSeg at the first time
    if si == 1
        idx_f = (fi>=f_AOI(1) &  fi<=f_AOI(2));       
        nf = length(find(idx_f));
        crossDensity_eachShuffle = zeros(shuffleN, nf);
        
        clear nf
    end
    phi1 = angle(S1(idx_f, :));
    phi2 = angle(S2(idx_f, :));
    
    
    crossDensity_eachShuffle(si, :) =  mean(exp(-1i * (phi1 - phi2)), 2);
    
    clear idx cross_density_sum
end

% abs imaginary of Coherency
psedo_iCoh = imag(crossDensity_eachShuffle);

% fit a normal distribution to psedo_iCoh
[~, nf] = size(psedo_iCoh);
mus = zeros(nf,1);
stds = zeros(nf, 1);
for fi = 1 : nf
    pd = fitdist(psedo_iCoh(:, fi),'Normal');
    mus(fi, 1) = pd.mu;
    stds(fi, 1) = pd.sigma;
    
    clear pd
    
end
