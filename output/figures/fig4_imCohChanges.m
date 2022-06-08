function fig4_imCohChanges()

codefilepath = mfilename('fullpath');


% find the codefolder
tmp = regexp(codefilepath, '.*\code', 'match');
if length(tmp) ~= 1
    disp('can not find code path correctly.')
    return;
end
codefolder = tmp{1};
clear tmp

% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));

%% Input & save
[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
[~, funcname, ~]= fileparts(codefilepath);


input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compCond');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compCond');


savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;



%% plot figure parameters
w_colormap = 350; % width  for the colormap
h_colormap = 120; % height for the colormap

w_deltax1_colormap = 5; % x distance between two color map within the same NHP
w_deltax2_colormap = 10; % x distance between two color map of different NHPs

w_textMovePhase = 70; % width showing the moveing phase, i.e. preMove
w_textpair = 80; % width showing the pair name, i.e. M1-STN
w_textColorbar = 80; % width showing the colarbar 

h_deltay_colormap = 15; % y distance between two color map 

h_textAnimal = 40; % height showing the animal, i.e. Animal J
h_textCond = 10; % height showing the condition, i.e. Mild-Normal
h_textFrenum = 30; % height showing the frequency number, i.e. 10 12
h_textFrelabel = 30; % height showing the frequency label, i.e. Frequences/Hz



fontsize1 = 11;
fontsize2 = 10;
fontname = 'Times New Roman';

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

fig_width = ncols * w_colormap + (ncols-2)* w_deltax1_colormap + w_deltax2_colormap + w_textpair + w_textMovePhase;
fig_height = nrows * h_colormap + (nrows-1)* h_deltay_colormap + h_textFrelabel + h_textFrenum + h_textCond + h_textAnimal;
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
    w_inner_left = 0;
    w_inner_right = 0;
    if coli == 1
        w_textpair_show = 0;
        w_inner_left = w_textpair;
        show_yticklabels = true;
    else
        w_textpair_show = w_textpair;
    end 
    if coli <= ncols_J
        w_outer_left = w_textMovePhase + w_textpair_show + (coli-1) * (w_colormap + w_deltax1_colormap);
    else
        w_outer_left = w_textMovePhase + w_textpair_show + (coli-1) * w_colormap + (coli-2) * w_deltax1_colormap + w_deltax2_colormap;
    end
    
    
    if coli == ncols
        w_textColorbar_show = 0;
        w_inner_right = w_textColorbar;
        show_colorbar = true;
    else
        w_textColorbar_show = w_textColorbar;
    end
    if coli > ncols_J
        w_outer_right = w_textColorbar_show + (ncols-coli) * (w_colormap + w_deltax1_colormap);
    else
        w_outer_right = w_textColorbar_show + (ncols-coli) * w_colormap + (ncols-coli -1) * w_deltax1_colormap + w_deltax2_colormap;
    end
    
    ifig_width = fig_width - (w_outer_left + w_outer_right);
    
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
        h_inner_top = 0;
        h_inner_bottom = 0;
        if rowi == 1
            h_textCond_show = 0;
            h_inner_top = h_textCond;
            show_titlename = true;
        else
            h_textCond_show = h_textCond;
        end
        h_outer_top = h_textAnimal + h_textCond_show + (rowi -1) * (h_colormap + h_deltay_colormap);
        if rowi == nrows || (rowi == nrows_Both && coli <= ncols_J)
            h_textFrenumlabel_show = 0;
            h_inner_bottom = h_textFrenum + h_textFrelabel;
            show_xlabel = true;
            show_xticklabels = true;
        else
            h_textFrenumlabel_show = h_textFrenum + h_textFrelabel;
        end
        h_outer_bottom = h_textFrenumlabel_show + (nrows-rowi)* (h_colormap + h_deltay_colormap);

        % outer and inner margin
        subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
        subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
        
        
        % actual plot
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, [comppd '-' basepd], 'histClim', [-1 1],...
            'codesavefolder', '', 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ...
            'fontsize1', fontsize1, 'fontsize2', fontsize2, 'fontname', fontname);
        
        % separate figure
        ifig_height = fig_height - (h_outer_top + h_outer_bottom);
        
        
        ifig = figure('Position', [50 50 ifig_width ifig_height]);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, [comppd '-' basepd], 'histClim', [-1 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig,  'innerposMargin', subplot_innerposMargin, ...
            'fontsize1', fontsize1, 'fontsize2', fontsize2, 'fontname', fontname);
        
        subfilename = [savefilename '-' animal '-' comppd '-B' basepd '-' event]; % 'Jo-mild-Bnormal-preMove'
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        close(ifig)
        
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

    t1 = annotation(fig, 'textbox', 'String', {event}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
    pos = t1.Position;
    pos_left = 5;
    pos_lower = fig_height-h_textAnimal-h_textCond-(rowi-1)*(h_colormap+h_deltay_colormap)-h_colormap/2-pos(4)*2/3;
    if rowi > nrows_Both
        pos_left = w_textMovePhase + w_textpair + w_colormap * ncols_J + w_deltax1_colormap * (ncols_J-1) + w_deltax2_colormap - pos(3)/2;
    end
    t1.Position = [pos_left pos_lower pos(3) pos(4)];
    
    clear event t pos lower left
end


%%%  added animal text
t1 = annotation(fig, 'textbox', 'String', {'Animal J'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
pos = t1.Position;
pos_left_J = (w_textMovePhase + w_textpair + w_colormap * ncols_J + w_deltax1_colormap * (ncols_J-1))/2;
pos_lower1 = fig_height-h_textAnimal - pos(4)/2;
t1.Position = [pos_left_J pos_lower1 pos(3) pos(4)];

t2 = annotation(fig, 'textbox', 'String', {'Animal K'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
pos = t2.Position;
pos_lower2 = pos_lower1;
pos_left_K = ((w_textMovePhase + w_textpair + w_colormap * ncols_J + w_deltax1_colormap * (ncols_J-1) + w_deltax2_colormap) + fig_width)/2 - pos(3)/2;
t2.Position = [pos_left_K pos_lower2 pos(3) pos(4)];


%%% save
print(fig, fullfile(savefolder, savefilename), '-painters', '-depsc');
print(fig, fullfile(savefolder, savefilename), '-dpng', '-r1000')
disp('saved figure')
close(fig)

