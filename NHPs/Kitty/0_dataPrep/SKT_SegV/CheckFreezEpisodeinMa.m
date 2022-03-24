load('H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\m3_fs500Hz_freezeSKTData_EpisodeExtract\Kitty_freezeEpisodes_moderate_20150408_bktdt2.mat')
savefolder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\m3_fs500Hz_freezeSKTData_EpisodeExtract';

freezEpisodes = freezStruct.freezEpisodes;

for frzi = 1: 5

    tri = freezEpisodes{frzi}.triali;
    freezeType = freezEpisodes{frzi}.freezeType;
    
    ma = smoothWspeed_trial{tri};
    tevents_ma = T_idxevent_ma{tri, :} / fs_ma;

    figure
    plot([1: length(ma)]/fs_ma, ma)
    hold on
    title(['tri = ' num2str(tri) ',' freezeType])

    for tei = 1 : length(tevents_ma)
        plot([tevents_ma(tei) tevents_ma(tei)], ylim, '--')
    end

    freeTphase = freezEpisodes{frzi}.freezeTPhaseS;
    for fti = 1 : length(freeTphase)
        ys = ylim;
        ys = [ys(1) ys(2)*2/3];
        plot([freeTphase(fti) freeTphase(fti)], ys)
    end
    saveas(gcf, fullfile(savefolder, [num2str(frzi) '-' freezeType '.tif']))
    close gcf
end
