function psedoFreeCiCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohPhasefile, freezType, varargin)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohPhasefile:
%       
%       Name-Value: 
%           'codesavefolder' - code saved folder         
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

nchns = size(lfptrials, 1);

load(ciCohPhasefile, 'psedociCohs');
if ~exist('psedociCohs', 'var')
    psedociCohs = struct();
end

if ~isfield(psedociCohs, freezType) 
    psedociCohs.(freezType) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs.(freezType), 4) + 1;
end
disp(['psedo test start at ' num2str(shuffi_str) ' times'])
for si = shuffi_str : suffi_end
    for chni = 1 : nchns -1
        lfp1 = squeeze(lfptrials(chni, :, :));
        for chnj = chni + 1 : nchns
            lfp2 = squeeze(lfptrials(chnj, :, :));
            [psedociCoh, ~] = psedo_ciCohSKT_FFT_NoAmp(lfp1, lfp2, fs, f_AOI);
            
            if chni == 1 && chnj == chni + 1
                nf = size(psedociCoh, 1);
                psedociCohs1Time = zeros(nchns, nchns, nf, 1);
                clear nf
            end
            
            psedociCohs1Time(chni, chnj, :, :) = psedociCoh;
            clear lfp2 psedociCoh
        end
        clear lfp1
    end
    psedociCohs.(freezType) = cat(4, psedociCohs.(freezType), psedociCohs1Time);
    clear psedociCohs1Time
    if(mod(si,10) == 0)
        disp(['psedo test finished ' num2str(si) ' times'])
        save(ciCohPhasefile, 'psedociCohs', '-append');
    end
end