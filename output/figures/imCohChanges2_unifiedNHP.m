w_colormap = 400; % width  for the colormap
h_colormap = 150; % height for the colormap
deltax1_colormap = 30; % x distance between two color map within the same NHP
deltax2_colormap = 50; % x distance between two color map of different NHPs
deltaxy_colormap = 20; % y distance between two color map 
w_textpair = 80; % width showing the pair name, i.e. M1-STN
w_textMovePhase = 80; % width showing the moveing phase, i.e. preMove
w_textColorbar = 80; % width showing the colarbar 
h_textFrenum = 10; % height showing the frequency number, i.e. 10 12
h_textFrelabel = 40; % height showing the frequency label, i.e. Frequences/Hz
h_textCond = 20; % height showing the condition, i.e. Mild-Normal
h_textAnimal = 30; % height showing the condition, i.e. Mild-Normal


nrows = 4;
ncols = 3;

fig_width = ncols * w_colormap + (ncols-2)* deltax1_colormap + deltax2_colormap + w_textpair + w_textMovePhase;
fig_height = nrows * h_colormap + (nrows-1)* deltaxy_colormap + h_textFrelabel + h_textFrenum + h_textCond + h_textAnimal;

fig = figure('Position', [50 50 fig_width fig_height]);
set(fig, 'PaperUnits', 'points');
for rowi = 1 : nrows
    show_titlename = false;
    show_xlabel = false;
    show_xticklabels = false;
    
    % extract outer_top, outer_bottom, inner_top and inner_bottom
    inner_top = 0;
    inner_bottom = 0;
    if rowi == 1
        h_textCond_show = 0;
        inner_top = h_textCond;
        show_titlename = true;
    else
        h_textCond_show = h_textCond;
    end
    outer_top = h_textAnimal + h_textCond_show + (rowi -1) * (h_colormap + deltaxy_colormap);
    if rowi == nrows
        h_textFrenumlabel_show = 0;
        inner_bottom = h_textFrenum + h_textFrelabel;
        show_xlabel = true;
        show_xticklabels = true;
    else
        h_textFrenumlabel_show = h_textFrenum + h_textFrelabel;
    end
    outer_bottom = h_textFrenumlabel_show + (nrows-rowi)* (h_colormap + deltaxy_colormap);
    
    for coli = 1 : ncols
        show_yticklabels = false;
        show_colorbar = false;
        
        % extract outer_left, outer_right, inner_left and inner_right
        inner_left = 0;
        inner_right = 0;
        if coli == 1
            w_textpair_show = 0;
            inner_left = w_textpair;
            show_yticklabels = true;
        else
            w_textpair_show = w_textpair;
        end
        outer_left = w_textMovePhase + w_textpair_show + (coli-1) * (w_colormap + deltax1_colormap);
        
        if coli == ncols
            w_textColorbar_show = 0;
            inner_right = w_textColorbar;
            show_colorbar = true;
        else
            w_textColorbar_show = w_textColorbar;
        end
        outer_right = w_textColorbar_show + (ncols-coli) * (w_colormap + deltax1_colormap);
        
        % outer and inner margin
        subplot_outerMargin = [outer_left outer_top outer_right outer_bottom];
        subplot_innerposMargin = [inner_left inner_top inner_right inner_bottom];
        
        %disp(['rowi = ' num2str(rowi) ', coli = ' num2str(coli)])
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, 'Mild-Normal', 'histClim', [-1 1],...
            'codesavefolder', '', 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin);
        
        clear subplot_outerMargin subplot_innerposMargin
        clear outer_left outer_right inner_left inner_right
        clear show_yticklabels show_colorbar
    end
    
    clear outer_top outer_bottom inner_top inner_bottom
    clear show_titlename show_xlabel show_xticklabels
end










