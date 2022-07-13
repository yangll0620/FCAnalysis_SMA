plot_fig1 = false;
plot_fig2 = true;
plot_fig3 = true;
plot_fig5 = false;



if plot_fig1
    plot_normalSpeed = true;
    plot_reachTime = true;
    
    if plot_normalSpeed
        fig1_normalTrial_SpeedExample('pos_ifig', [150 150 500 120]);
    end

    if plot_reachTime
        fig1_reachTimeStatiscal('plot_timeStatiscal', true, 'pos_ifig_reachTimeStatiscal', [150 150 300 200])
    end
end

if plot_fig2
    fig2_imCohChanges_compCond('pos_ifig', [50 50 400 200])
end

if plot_fig3
    fig3_imCohChanges_compEvent('pos_ifig', [50 50 400 200])
end


if plot_fig5
    plot_freezeSpeed = true;
    plot_freezeImCoh = false;
    
    if plot_freezeSpeed
        fig5_freeze_SpeedExample('pos_ifig', [150 150 400 200])
    end

    if plot_freezeImCoh
        fig5_freeze_imCoh('pos_ifig', [50 50 400 200])
    end
end