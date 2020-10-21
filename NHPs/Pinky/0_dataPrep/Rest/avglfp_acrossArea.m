function [avglfp, T_chnsarea_new] = avglfp_acrossArea(lfpdata, T_chnsarea)
    % average lfp across each area 
    % 
    % args:
    %       lfpdata : lfp data in area (nchns * ntemp)
    %       T_chnsarea: chns-area table 
    % return:
    %       avglfp: averaged lfp data (nareas * ntemp)
    %       T_chnsarea_new: the new chns-area table


    uniqAreas = unique(T_chnsarea.brainarea);
    avglfp = [];
    T_chnsarea_new = T_chnsarea([], :);
    for ai = 1 : length(uniqAreas)
        barea = uniqAreas{ai};


        mask_area = strcmp(T_chnsarea.brainarea, barea);
        lfp_area = lfpdata(mask_area, :);
        avglfp_area = mean(lfp_area, 1);


        avglfp = cat(1, avglfp, avglfp_area);

        tbl = T_chnsarea(find(mask_area,1),:);
        tbl.recordingchn = nan;
        T_chnsarea_new = [T_chnsarea_new; tbl];

        clear barea mask_area avglfp_area tbl
    end
    T_chnsarea_new.chni = [1: height(T_chnsarea_new)]';
end