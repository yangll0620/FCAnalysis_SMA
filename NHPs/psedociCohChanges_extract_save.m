function psedociCohChanges_extract_save(suffi_end, lfptrials, lfpsegs_base, fs, f_AOI, ciCohChangesfile)
%
% lfptrials: nchns * ntemp * ntrials(nsegs)
%   
% save psedoiCohChanges: nchns * nchns * nf * nshuffle, saved to ciCohChangesfile

load(ciCohChangesfile, 'psedoiCohChanges');

if(~exist('psedoiCohChanges', 'var'))
    psedoiCohChanges = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedoiCohChanges, 4) + 1;
end
lfp_combined = cat(3, lfpsegs_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_SKT = lfp_combined(:, :, masksBase);
    psedolfp_Base = lfp_combined(:, :, ~masksBase);
    
    [~, psedoiCoh_SKT, ~] = ciCoh_trialDeltaPhi(psedolfp_SKT, fs, f_AOI);
    [~, psedoiCoh_Base, ~] = ciCoh_trialDeltaPhi(psedolfp_Base, fs, f_AOI);
    
    psedoiCohChanges = cat(4, psedoiCohChanges, psedoiCoh_SKT - psedoiCoh_Base);
    
    if(mod(si, 10) == 0)
        disp(['pesdo test finished ' num2str(si)])
        save(ciCohChangesfile, 'psedoiCohChanges', '-append');
    end
    
    clear randomSKTInds randomRestInds psedolfp_SKT psedoiCoh_rest
    clear psedoiCoh_SKT psedoiCoh_rest
end