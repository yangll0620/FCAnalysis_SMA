function [eqLen_avglfp, T_chnsarea_new, fs] = combineVLVPL_1file(file, segt)
    %%   segment into intervals with same length and average lfp across each area for data in file
    %       combine VLo and VPLo
    %
    % Args:
    %   file:  one file (full path)
    %   segt: the same time length, unit second (e.g 2)
    %
    %  Returns:
    %   eqLen_avglfp: averaged lfp with same length ((nareas + nDBS) * ntemp * nsegs)
    %   T_chnsarea_new: new T_chnsarea_new (height = (nareas + nDBS))
    %   fs: sample rate

    load(file, 'data_segments', 'fs', 'T_chnsarea')
    
   
     
end