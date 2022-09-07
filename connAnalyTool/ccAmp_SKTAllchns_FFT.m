function [ccAmp, f_selected] = ccAmp_SKTAllchns_FFT(lfpdata, fs, f_AOI, varargin)
%
%   ref: https://yll0620.medium.com/functional-connectivity-measurement-16423fee3581 
%      treat ntemp as stationary time series, i.e. twin = ntemp, toverlap = 0
%
%   Input
%       lfpdata: nchns * ntemp * ntrials
%       
%       fs: sampling rate
%
%       f_AOI: frequencies duration of interest, i.e. f_AOI = [8 40]
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
%
%   Outputs:
%       ccAmp: correlation coefficient between amplitudes (nchns * nchns * nf), Symmetrical



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

[nchns, ~, ntrials] = size(lfpdata);

% extract amps_allchns: nchns * nf * ntrials
for chni = 1 : nchns
    lfptrials = squeeze(lfpdata(chni, :, :));
    [~, amps, f_selected] = phaseAmp_SKTPerTrial_FFT(lfptrials, fs, f_AOI);
    
    if chni == 1
        [nf, ~] = size(amps);
        amps_allchns = zeros(nchns, nf, ntrials); 
        clear nf 
    end
    
    amps_allchns(chni, :, :) = amps;
    
    clear lfptrials amps  
end


% calculate ccAmp
nf = size(amps_allchns, 2);
ccAmp = zeros(nchns, nchns, nf);
for nfi = 1 : nf
    for chni = 1: nchns -1
        ampi = squeeze(amps_allchns(chni, nfi, :));
        
        for chnj = chni + 1 : nchns
            ampj = squeeze(amps_allchns(chnj, nfi, :));
            
            [r, p] = corrcoef(ampi,ampj);
            if(p(1,2) >= 0.05) % not sig (P is smaller than the significance level (default is 0.05), then the corresponding correlation in R is considered significant.)
                r(1,2) = 0;
            end
            ccAmp(chni, chnj, nfi) = r(1,2);
        end
        clear ampi
    end
end
