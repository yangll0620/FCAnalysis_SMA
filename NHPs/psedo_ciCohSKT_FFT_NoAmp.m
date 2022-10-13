function [psedociCoh_NoAmp, f_selected] = psedo_ciCohSKT_FFT_NoAmp(lfp1, lfp2, fs, f_AOI, varargin)
%
%   ref: https://yll0620.medium.com/functional-connectivity-measurement-16423fee3581 
%      treat ntemp as stationary time series, i.e. twin = ntemp, toverlap = 0
%
%   Input
%       lfp1, lfp2: ntemp * ntrials
%
%       f_AOI: frequencies duration of interest, i.e. f_AOI = [8 40]
%
%       Name-Value: 
%           'codesavefolder' - code saved folder 
%
% return:
%       psedociCoh_NoAmp: psedo absolute imaginery coherence nf * 1


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

% shuffle lfp2
ntrials = size(lfp2, 2);
lfp2_psedo = lfp2(:, randperm(ntrials));


[phis1, ~, f_selected] = phaseAmp_SKTPerTrial_FFT(lfp1, fs, f_AOI);
[phis2, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfp2_psedo, fs, f_AOI);

deltaphi = phis1 - phis2;
coh = mean(exp(1i * deltaphi), 2);
psedociCoh_NoAmp = abs(imag(coh) ./ sqrt((1 - real(coh).^2)));