clear 

%%function imCohChanges_unifiedNHP()



% add path
codefolder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\code';
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));

%% Input
input_folder_J = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Jo\0_dataPrep\SKT\fs500Hz\m4_fs500Hz_uNHP_imCohChanges_compCond';
input_folder_K = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT\fs500Hz\m4_fs500Hz_uNHP_imCohChanges_compCond';
savecodefolder = '';


%% plot figure parameters
w_colormap = 400; % width  for the colormap
h_colormap = 120; % height for the colormap
deltax1_colormap = 10; % x distance between two color map within the same NHP
deltax2_colormap = 15; % x distance between two color map of different NHPs
deltaxy_colormap = 20; % y distance between two color map 
w_textpair = 80; % width showing the pair name, i.e. M1-STN
w_textMovePhase = 60; % width showing the moveing phase, i.e. preMove
w_textColorbar = 80; % width showing the colarbar 
h_textFrenum = 10; % height showing the frequency number, i.e. 10 12
h_textFrelabel = 40; % height showing the frequency label, i.e. Frequences/Hz
h_textCond = 20; % height showing the condition, i.e. Mild-Normal
h_textAnimal = 30; % height showing the condition, i.e. Mild-Normal


%% Code start here
basepd = 'normal';

conds_J = cond_cell_extract('Jo');
conds_K = cond_cell_extract('Kitty');
conds_J(strcmp(conds_J, basepd)) = [];
conds_K(strcmp(conds_K, basepd)) = [];
conds = [conds_J conds_K];


ePhases_both = {'preMove'; 'earlyReach';  'lateReach'};
ePhases_onlyK = {'PeakV'};
ePhases = [ePhases_both; ePhases_onlyK];

nrows = length(ePhases);
nrows_Both = length(ePhases_both);
ncols_J = length(conds_J);
ncols = length(conds);

fig_width = ncols * w_colormap + (ncols-2)* deltax1_colormap + deltax2_colormap + w_textpair + w_textMovePhase;
fig_height = nrows * h_colormap + (nrows-1)* deltaxy_colormap + h_textFrelabel + h_textFrenum + h_textCond + h_textAnimal;
fig = figure('Position', [50 50 fig_width fig_height]);
set(fig, 'PaperUnits', 'points');
for coli = 1 : ncols
    % extract comppd and animal
    comppd = conds{coli};
    if coli <= ncols_J
        animal = 'Jo';
        input_folder = input_folder_J;
    else
        animal = 'Kitty';
        input_folder = input_folder_K;
    end
    ciCohChangesfile_prefix =[animal '-ciCohChanges'];
    
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
    if coli <= ncols_J
        outer_left = w_textMovePhase + w_textpair_show + (coli-1) * (w_colormap + deltax1_colormap);
    else
        outer_left = w_textMovePhase + w_textpair_show + (coli-1) * w_colormap + (coli-2) * deltax1_colormap + deltax2_colormap;
    end
    
    
    if coli == ncols
        w_textColorbar_show = 0;
        inner_right = w_textColorbar;
        show_colorbar = true;
    else
        w_textColorbar_show = w_textColorbar;
    end
    if coli > ncols_J
        outer_right = w_textColorbar_show + (ncols-coli) * (w_colormap + deltax1_colormap);
    else
        outer_right = w_textColorbar_show + (ncols-coli) * w_colormap + (ncols-coli -1) * deltax1_colormap + deltax2_colormap;
    end
    
    
    for rowi = 1 : nrows
        
        % extract sigciCohChanges_flatten
        event = ePhases{rowi};
        [~, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, comppd, 'codesavefolder', savecodefolder);
        ciCohChangesfile = fullfile(input_folder, [ciCohChangesfile_prefix  '_b' basepd '--' comppd '_' event '_align2' align2name '.mat']);
        if ~exist(ciCohChangesfile, 'file')
            clear event align2name ciCohChangesfile
            continue;
        end
        load(ciCohChangesfile, 'ciCohChanges', 'psedoiCohChanges', 'f_selected',  'T_chnsarea')
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
    
        
        
        % extract outer_top, outer_bottom, inner_top and inner_bottom and set show tag
        show_titlename = false;
        show_xlabel = false;
        show_xticklabels = false;
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

        % outer and inner margin
        subplot_outerMargin = [outer_left outer_top outer_right outer_bottom];
        subplot_innerposMargin = [inner_left inner_top inner_right inner_bottom];
        
        
        % actual plot
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, [comppd '-' basepd], 'histClim', [-1 1],...
            'codesavefolder', '', 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin);
        
        clear subplot_outerMargin subplot_innerposMargin
        clear outer_top outer_bottom inner_top inner_bottom
        clear show_titlename show_xlabel show_xticklabels
        clear event align2name ciCohChangesfile
        clear ciCohChanges psedoiCohChanges f_selected  T_chnsarea
        clear sigciCohChanges sigciCohChanges_flatten chnPairNames
    end
    
    clear outer_left outer_right inner_left inner_right
    clear show_yticklabels show_colorbar
    clear comppd ciCohChangesfile_prefix
    
end


%%%  added event text
for rowi = 1 : nrows
    event = ePhases{rowi};

    t = annotation(fig, 'textbox', 'String', {event}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 11, 'FontWeight', 'bold');
    pos = t.Position;
    pos_left = 5;
    pos_lower = fig_height-h_textAnimal-h_textCond-(rowi-1)*(h_colormap+deltaxy_colormap)-h_colormap/2-pos(4)*2/3;
    if rowi > nrows_Both
        pos_left = w_textMovePhase + w_textpair + w_colormap * ncols_J + deltax1_colormap * (ncols_J-1) + deltax2_colormap - pos(3)/2;
    end
    t.Position = [pos_left pos_lower pos(3) pos(4)];
    
    clear event t pos lower left
end


%%%  added animal text
pos_left_J = (w_textMovePhase + w_textpair + w_colormap * ncols_J + deltax1_colormap * (ncols_J-1))/2;
pos_left_K = ((w_textMovePhase + w_textpair + w_colormap * ncols_J + deltax1_colormap * (ncols_J-1) + deltax2_colormap) + fig_width)/2; 
pos_lower = fig_height-h_textAnimal;

t = annotation(fig, 'textbox', 'String', {'Animal J'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 11, 'FontWeight', 'bold');
pos = t.Position;
t.Position = [pos_left_J pos_lower pos(3) pos(4)];

t = annotation(fig, 'textbox', 'String', {'Animal K'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 11, 'FontWeight', 'bold');
pos = t.Position;
t.Position = [pos_left_K pos_lower pos(3) pos(4)];
