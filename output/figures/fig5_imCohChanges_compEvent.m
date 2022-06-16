function fig5_imCohChanges_compEvent()
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

input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');


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

h_deltay_colormap_J = 80; % y distance between two color map of animal J
h_deltay_colormap_K = 10; % y distance between two color map of animal K

h_textAnimal = 40; % height showing the animal name, i.e. animal J/K
h_textCond = 10; % height showing the condition, i.e. Mild-Normal
h_textFreNum = 30; % height showing the frequency number, i.e. 10 12
h_textFreLabel = 30; % height showing the frequency label, i.e. Frequences/Hz


fontname = 'Times New Roman';


%% Code start here
baseevent = 'preMove';

ePhases_J = {'earlyReach';  'lateReach'};
ePhases_K = {'earlyReach';  'PeakV'; 'lateReach'};

cond_J = {'Normal'; 'Mild';  'Moderate'};
cond_K = {'Normal'; 'Moderate'};

nrows_J = length(ePhases_J);
nrows_K = length(ePhases_K);
ncols_J = length(cond_J);
ncols_K = length(cond_K);

ncols = ncols_J + ncols_K;

close all
fig_width = ncols * w_colormap + (ncols-2)* w_deltax1_colormap + w_deltax2_colormap + w_textMovePhase * 2 + w_textpair * 2 + w_textColorbar * 2;
fig_height_J = h_textAnimal + + h_textCond + nrows_J * h_colormap + (nrows_J-1)* h_deltay_colormap_J + h_textFreLabel + h_textFreNum;
fig_height_K = h_textAnimal + + h_textCond + nrows_K * h_colormap + (nrows_K-1)* h_deltay_colormap_K + h_textFreLabel + h_textFreNum;
fig_height = max(fig_height_J, fig_height_K);
fig = figure('Position', [50 50 fig_width fig_height]);
set(fig, 'PaperUnits', 'points');


%%%  added animal text
t1 = annotation(fig, 'textbox', 'String', {'Animal J'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
pos = t1.Position;
pos_left_J = (w_textMovePhase + w_textpair + w_colormap * ncols_J + w_deltax1_colormap * (ncols_J-1))/2;
pos_lower1 = fig_height-h_textAnimal - pos(4)/2;
t1.Position = [pos_left_J pos_lower1 pos(3) pos(4)];

t2 = annotation(fig, 'textbox', 'String', {'Animal K'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
pos = t2.Position;
pos_lower2 = pos_lower1;
pos_left_K = ((w_textMovePhase*2 + w_textpair*2 + w_colormap * ncols_J + w_deltax1_colormap * (ncols_J-1) + w_deltax2_colormap + w_textColorbar) + fig_width)/2 - pos(3)/2;
t2.Position = [pos_left_K pos_lower2 pos(3) pos(4)];


%%%  added event text
for rowi = 1 : nrows_J
    event = ePhases_J{rowi};

    t1 = annotation(fig, 'textbox', 'String', {event}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
    pos = t1.Position;
    pos_left = 5;
    pos_lower = fig_height-h_textAnimal-h_textCond-(rowi-1)*(h_colormap+h_deltay_colormap_J)-h_colormap/2-pos(4)*2/3;
    t1.Position = [pos_left pos_lower pos(3) pos(4)];
    
    clear event t pos lower left
end
for rowi = 1 : nrows_K
    event = ePhases_K{rowi};

    t1 = annotation(fig, 'textbox', 'String', {event}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
    pos = t1.Position;
    pos_left = w_textMovePhase + w_textpair + ncols_J * w_colormap  + (ncols_J-1) * w_deltax1_colormap + w_deltax2_colormap + w_textColorbar;
    pos_lower = fig_height-h_textAnimal-h_textCond-(rowi-1)*(h_colormap+h_deltay_colormap_K)-h_colormap/2-pos(4)*2/3;
    t1.Position = [pos_left pos_lower pos(3) pos(4)];
    
    clear event t pos lower left
end

%%% plot actual changes hist plot
for coli = 1 : ncols
    
    if coli <= ncols_J % for animal J
        animal = 'Jo';
        input_folder = input_folder_J;
        pdcond = cond_J{coli};
        nrows = nrows_J;
        ePhases = ePhases_J;
        h_deltay_colormap = h_deltay_colormap_J;
    else
        animal = 'Kitty';
        input_folder = input_folder_K;
        pdcond = cond_K{coli-ncols_J};
        nrows = nrows_K;
        ePhases = ePhases_K;
        h_deltay_colormap = h_deltay_colormap_K;
    end

    show_textpair = false;
    show_colorbar = false;
    if coli == 1 || coli == ncols_J +1
        show_textpair = true;
    end
    if coli == ncols_J || coli == ncols
        show_colorbar = true;
    end

    if coli <= ncols_J
        w_outer_left = w_textMovePhase + w_textpair + (coli - 1) * (w_colormap + w_deltax1_colormap);
    else
        w_outer_left = w_textMovePhase * 2 + w_textpair * 2 + (coli - 1) * w_colormap + (coli-2) * w_deltax1_colormap + w_textColorbar + w_deltax2_colormap;
    end
    w_outer_diff = w_colormap;
    w_inner_left = 0;
    w_inner_right = 0;
    if show_textpair
        w_outer_left = w_outer_left - w_textpair;
        w_outer_diff = w_outer_diff + w_textpair;
        w_inner_left = w_inner_left + w_textpair;
    end
    if show_colorbar
        w_outer_diff = w_outer_diff + w_textColorbar;
        w_inner_right = w_inner_right + w_textColorbar;
    end
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);
    
    for rowi = 1 : nrows
        
        % extract sigciCohChanges_flatten
        event = ePhases{rowi};
        [~, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, 'codesavefolder', savecodefolder);
        pattfile = dir(fullfile(input_folder, [animal  '*ciCohChanges*' lower(pdcond) '*b' baseevent '--' event '_align2' align2name '.mat']));
        if length(pattfile) ~=1
            disp('exist file is not only one')
            return;
        end
        ciCohChangesfile = fullfile(input_folder, pattfile(1).name);
        load(ciCohChangesfile, 'ciCohChanges', 'psedoiCohChanges', 'f_selected',  'T_chnsarea')
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        
        
        % extract outer_top, outer_bottom, inner_top and inner_bottom and set show tag
        show_condname = false;
        show_freLabel = false;
        show_freNum = false;
        
        if rowi == 1
            show_condname = true;
        end
        if rowi == nrows
            show_freLabel = true;
            show_freNum = true;
        end

        
        h_outer_top = h_textAnimal + h_textCond + (rowi - 1) * (h_colormap + h_deltay_colormap);
        h_outer_diff = h_colormap;
        h_inner_top = 0;
        h_inner_bottom = 0;
        if show_condname
            h_outer_diff = h_outer_diff + h_textCond;
            h_outer_top = h_outer_top - h_textCond;
            h_inner_top = h_inner_top + h_textCond;
        end
        if show_freLabel
            h_outer_diff = h_outer_diff + h_textFreLabel;
            h_inner_bottom =  h_inner_bottom + h_textFreLabel;
        end
        if show_freNum
            h_outer_diff = h_outer_diff + h_textFreNum;
            h_inner_bottom =  h_inner_bottom + h_textFreNum; 
        end
        h_outer_bottom = fig_height - (h_outer_top + h_outer_diff);

        
        % outer and inner margin
        subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
        subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
        
        % actual plot
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, pdcond, 'histClim', [-1 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_freNum, 'show_yticklabels', show_textpair, 'show_xlabel', show_freLabel, 'show_titlename', show_condname,'show_colorbar', show_colorbar, ...
            'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ...
            'fontname', fontname);
        
        ifig = figure('Position', [150 150 400 200]);
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, pdcond, 'histClim', [-1 1],...
            'fig', ifig, ...
            'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_freNum, 'show_yticklabels', show_textpair, 'show_xlabel', show_freLabel, 'show_titlename', show_condname,'show_colorbar', show_colorbar, ...
            'codesavefolder', savecodefolder);
        subfilename = [animal '-' event '-' pdcond]; % 'Jo-mild-Bnormal-preMove'
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
        close(ifig)
        
    end
end

%%% save
print(fig, fullfile(savefolder, savefilename), '-painters', '-depsc');
print(fig, fullfile(savefolder, savefilename), '-dpng', '-r1000')
disp('saved figure')
close(fig)

