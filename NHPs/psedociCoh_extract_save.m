function psedociCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohPhasefile)
%
% lfptrials: nchns * ntemp * ntrials(nsegs)
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile

nchns = size(lfptrials, 1);

load(ciCohPhasefile, 'psedociCohs');
if(~exist('psedociCohs', 'var'))
    psedociCohs = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs, 4) + 1;
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
    psedociCohs = cat(4, psedociCohs, psedociCohs1Time);
    clear psedociCohs1Time
    if(mod(si,10) == 0)
        disp(['psedo test finished ' num2str(si) ' times'])
        save(ciCohPhasefile, 'psedociCohs', '-append');
    end
end