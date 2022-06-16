function [phis, amps, f_selected] = phaseAmp_FFT_extract(lfptrials, fs, f_AOI, varargin)
% 
%   Extract phase and Amplitudes using FFT
%
%   Input
%       lfpdata:  ntemp * ntrials
%
%       f_AOI: frequencies duration of interest, i.e. f_AOI = [8 40]
%
%   Outputs:
%       phis: phase angle in radians for each trial (nf * ntrials)
%
%       amps: magnitude for each trial (nf * ntrials)
%
%       f_selected: selected frequencies (nf * 1)


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

codesavefolder = p.Results.codesavefolder;

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


[ntemp, ntrials] = size(lfptrials);
phis = [];
amps = [];
for triali = 1: ntrials
    x = lfptrials(: , triali);
    
    [Sx, fx, ~, ~] = spectrogram(x, ntemp, 0,[],fs); % Sx: nf * 1
    
    if triali == 1
        freqs = fx;
        idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
        f_selected = round(freqs(idx_f), 2);
        
        clear freqs
    end
    
    phi = angle(Sx(idx_f, :)); % phix : nf *1
    amp = abs(Sx(idx_f, :)); % amp: nf * 1
    
    phis = cat(2, phis, phi);
    amps = cat(2, amps, amp);
    
    clear x Sx fx phi amp
end
