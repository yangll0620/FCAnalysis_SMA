function psedociCohChanges_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile, varargin)
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
% save psedoiCohChanges: nchns * nchns * nf * nshuffle, saved to ciCohChangesfile


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


load(ciCohChangesfile, 'psedoiCohChanges');

if(~exist('psedoiCohChanges', 'var'))
    psedoiCohChanges = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedoiCohChanges, 4) + 1;
end
lfp_combined = cat(3, lfptrials_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_SKT = lfp_combined(:, :, masksBase);
    psedolfp_Base = lfp_combined(:, :, ~masksBase);
    
    [~, psedoiCoh_SKT, ~] = ciCoh_trialDeltaPhi(psedolfp_SKT, fs, f_AOI, 'codesavefolder', codesavefolder);
    [~, psedoiCoh_Base, ~] = ciCoh_trialDeltaPhi(psedolfp_Base, fs, f_AOI);
    
    psedoiCohChanges = cat(4, psedoiCohChanges, psedoiCoh_SKT - psedoiCoh_Base);
    
    if(mod(si, 10) == 0)
        disp(['pesdo test finished ' num2str(si)])
        save(ciCohChangesfile, 'psedoiCohChanges', '-append');
    end
    
    clear randomSKTInds randomRestInds psedolfp_SKT psedoiCoh_rest
    clear psedoiCoh_SKT psedoiCoh_rest
end