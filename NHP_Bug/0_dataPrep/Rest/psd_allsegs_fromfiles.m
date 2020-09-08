function [psd_allsegs, F_pxx] = psd_allsegs_fromfiles(files, twin_pwelch, brainarea)
%% extract psd of all the segments from all files
%
% Outputs:
%   
%       psd_allsegs: the PSD estimate of all the segments from the files  (nfs * nsegs)
%
%       F_pxx: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
%
%       brainarea: the brain area to analysis ('stn', or 'gp')




nfiles = length(files);
psd_allsegs = [];
for filei = 1 : nfiles
    
    % load data
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'data_segments', 'fs', 'segsRemain');
    
    
    for segi = 1: length(segsRemain)
        if segsRemain(segi) == 0 % ignore the segment marked with 0
            continue;
        end
        
        
        % extract the lfp data of segi in brainarea
        eval(['lfp_oneseg = data_segments(segi).lfp_' brainarea ';']);
        
        
        %%% calcualte PSD throuch zscore and pwelch%%%
        
        % zscore of the lfp_oneseg along each column
        lfp_zscore = zscore(lfp_oneseg);
        
        
        % psd pwelch paramers
        nwins = round(twin_pwelch * fs);
        noverlap = round(nwins * 0.9);
        
        % PSD is computed independently for each column, Pxx: nfs * nchns, f: nfs * 1
        [Pxx, f]= pwelch(lfp_zscore, nwins, noverlap, nwins, fs);
        
        
        if ~exist('F_pxx', 'var')
            F_pxx = f;
        else
            if ~isequal(F_pxx, f)
                disp(['F_pxx not equal f for ' filename ', segi = ' num2str(segi)])
                continue;
            end
        end
        
        % psd_allsegs: nfs * nchns * nsegs
        psd_allsegs = cat(3, psd_allsegs, Pxx);
        
        clear lfp_oneseg lfp_zscore Pxx f
    end
    
    clear filename data_segments fs segsRemain
end