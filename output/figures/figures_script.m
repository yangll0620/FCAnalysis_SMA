plot_fig1 = false;
plot_fig5 = true;



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


if plot_fig5
    plot_freezeSpeed = true;
    if plot_freezeSpeed
        fig5_freeze_SpeedExample('pos_ifig', [150 150 500 250])
    end
end