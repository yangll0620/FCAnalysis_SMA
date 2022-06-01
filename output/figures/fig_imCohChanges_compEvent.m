function fig_imCohChanges_compEvent()
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
input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');


savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;


%% plot figure parameters
w_colormap = 350; % width  for the colormap
h_colormap = 120; % height for the colormap

w_deltax1_colormap = 5; % x distance between two color map within the same NHP
w_deltax2_colormap = 20; % x distance between two color map of different NHPs

w_textMovePhase = 70; % width showing the moveing phase, i.e. preMove
w_textpair = 80; % width showing the pair name, i.e. M1-STN
w_textColorbar = 80; % width showing the colarbar 

h_deltay_colormap_J = 80; % y distance between two color map of animal J
h_deltay_colormap_K = 10; % y distance between two color map of animal K

h_textAnimal = 30; % height showing the animal name, i.e. animal J/K
h_textCond = 30; % height showing the condition, i.e. Mild-Normal
h_textFreNum = 10; % height showing the frequency number, i.e. 10 12
h_textFreLabel = 40; % height showing the frequency label, i.e. Frequences/Hz


fontname = 'Times New Roman';


%% Code start here
baseevent = 'preMove';

ePhases_J = {'earlyReach';  'lateReach'};
ePhases_K = {'earlyReach';  'PeakV'; 'lateReach'};

cond_J = {'Mild';  'Moderate'};
cond_K = {'Moderate'};

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
        w_outer_left = w_textMovePhase * 2 + w_textpair * 2 + (coli - 1) * w_colormap + (ncols_J-2) * w_deltax1_colormap + w_textColorbar + w_deltax2_colormap;
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
        plot_ciCohHistogram3(sigciCohChanges_flatten, chnPairNames, f_selected, pdcond, 'histClim', [-1 1],...
            'codesavefolder', '', 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_freNum, 'show_yticklabels', show_textpair, 'show_xlabel', show_freLabel, 'show_titlename', show_condname,'show_colorbar', show_colorbar, ...
            'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ...
            'fontname', fontname);
        
    end
end

%%% save
print(fullfile(savefolder, savefilename), '-dpng', '-r1000')
disp('saved figure')


function plot_ciCohHistogram3(ciCoh_flatten, chnPairNames, f_selected, titlename, varargin)
%
%   Inputs:
%       ciCoh_flatten:
%       chnPairNames
%       f_selected
%       titlename
%       histClim
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
%           'histClim' - ciCoh histogram clim (default [0 1])
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'cbarTicks' - vector, default [0 0.5 1]
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_titlename' - show (true, default) or not show (false) titlename
%           'show_colorbar' - show (true, default) or not show (false) colorbar
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default [5 5 5 5]
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default [5 5 5 5]
%           'fontsize1' - font size for title, default 12
%           'fontsize2' - font size for , default 10
%           'fontname' - font name , default 'Times New Roman'



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [0 1], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'cbarTicks', [0 0.5 1], @(x) assert(isvector(x) && isnumeric(x)));
addParameter(p, 'cbarStr', 'ciCoh', @isstr);
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_titlename', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'innerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'fontsize1', 11, @(x) assert(isscalar(x) && isnumeric(x)));
addParameter(p, 'fontsize2', 10, @(x) assert(isscalar(x) && isnumeric(x)));
addParameter(p, 'fontname', 'Times New Roman', @isstr);

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig = p.Results.fig;
cbarTicks = p.Results.cbarTicks;
cbarStr = p.Results.cbarStr;
show_xticklabels = p.Results.show_xticklabels;
show_yticklabels = p.Results.show_yticklabels;
show_xlabel = p.Results.show_xlabel;
show_titlename = p.Results.show_titlename;
show_colorbar = p.Results.show_colorbar;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin; 
fontsize1 = p.Results.fontsize1;
fontsize2 = p.Results.fontsize2;
fontname = p.Results.fontname;



% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% plot
if isempty(fig)
    fig = figure;
end
ax = axes(fig, 'Units', 'pixels');
imagesc(ax, ciCoh_flatten)
colormap(jet)

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
set(gca,'CLim', histClim)


%%% plot the line to separete the pairs
chnPair_prev = '';
for ci = 1: length(chnPairNames)
    chnPair = chnPairNames{ci};
    
    % replace M1-stn0-1 to M1-STN
    s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
    if ~isempty(s_stn)
        for si = 1 : length(s_stn)
            chnPair = strrep(chnPair, s_stn{si}, 'STN');
        end
    end
    % replace M1-stn0-1 to M1-STN
    s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
    if ~isempty(s_gp)
        for si = 1 : length(s_gp)
            chnPair = strrep(chnPair, s_gp{si}, 'GP');
        end
    end
    
    if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
        hold on; plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'w--')
        % Create line
    end
    chnPair_prev = chnPair;
    chnPairNames{ci} = chnPair;
    
    
    clear s_stn s_gp chnPair
end



%%% show inf
[npairs, nf] = size(ciCoh_flatten);
if show_xticklabels
    
    xticks([1:nf])
    xticklabels(round(f_selected))
    
    set(gca, 'fontsize',fontsize2, 'FontName', fontname)
else
    xticks([]);
end

if show_yticklabels
    yticks([1:npairs]);
    set(gca,'YTickLabel',chnPairNames,'fontsize',11, 'FontName', fontname)
else 
    yticks([]);
end

if show_xlabel
    xlabel('Frequency (Hz)', 'fontsize', 12, 'FontName', fontname)
end
if show_titlename
    title(titlename, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_colorbar
    pos = get(gca, 'Position');
    c = colorbar;
    c.Label.String = cbarStr;
    if ~isempty(cbarTicks)
        set(c, 'Ticks', cbarTicks, 'FontSize', fontsize1, 'FontWeight', 'bold', 'FontName', fontname)
    end
    set(gca, 'Position', pos);
    clear pos
end


