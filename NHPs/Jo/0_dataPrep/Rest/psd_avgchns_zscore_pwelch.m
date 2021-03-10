function [Pxx,F]  = psd_avgchns_zscore_pwelch(lfpdata, twin, fs)
%% using pwelch to calculate power spectral density (PSD)
%
%   Steps:
%       1. average across channels
%       2. calulate the zscore of the averaged lfp
%       3. psd estiamtes using pwelch with twin, toverlap = twin * 0.9 and fs 
%
%   
%
%   Inputs:
%
%       lfpdata (ntemp * nchns) : lfp data 
%       
%       twin: the time length (s) for pwelch
%
%  Outputs:
%   
%       Pxx: the PSD estimate   (nfs * 1)
%
%       F: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)


%% Codes start here

% average across channels
lfp_avg = mean(lfpdata,2);


% zscore of the averaged lfp
lfp_zscore = zscore(lfp_avg);


% psd estimates via pwelch
nwins = round(twin * fs);
noverlap = round(nwins * 0.9);

[Pxx,F] = pwelch(lfp_zscore, nwins, noverlap, nwins, fs);
