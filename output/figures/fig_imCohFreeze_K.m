function fig_imCohFreeze_K()
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
input_folder = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_FreezeSKT_earlyMiddLateFreeze_imCohUsingFFT');


savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;


%% plot figure parameters
w_colormap = 350; % width  for the colormap
h_colormap = 120; % height for the colormap


w_textpair = 80; % width showing the pair name, i.e. M1-STN
w_textColorbar = 80; % width showing the colarbar 

h_deltay_colormap = 10; % y distance between two color map of animal K

h_textFreeze = 30; % height showing the animal name, i.e. animal J/K
h_textFreNum = 10; % height showing the frequency number, i.e. 10 12
h_textFreLabel = 40; % height showing the frequency label, i.e. Frequences/Hz


%% Code start here
freezeTypes = {'earlyFreeze', 'middleFreeze', 'lateFreeze'};
nfTypes = length(freezeTypes);
fig_height = h_colormap * nfTypes + h_deltay_colormap * (nfTypes -1) + h_textFreeze * nfTypes + (h_textFreNum + h_textFreLabel);
fig_width = w_textpair + w_colormap + w_textColorbar;
fig = figure('Position', [50 50 fig_width fig_height]);
set(fig, 'PaperUnits', 'points');

for rowi = 1 : nfTypes
    freeType = freezeTypes {rowi};
    
    ciCohPhasefile = ['Kitty-Freeze-' freeType '_ciCoh_ntrial84_moderate.mat'];
    load(fullfile(input_folder, ciCohPhasefile), 'ciCohs', 'T_chnsarea', 'nsegs', 'f_selected', 'psedociCohs');
    
    [sigciCoh]= sigciCoh_extract(psedociCohs.combinedfreezTypes, ciCohs.combinedfreezTypes);
    [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);
    
    

    show_freqNumLabel= false;
    show_pairname = true;
    show_freezeType = true;
    if rowi == nfTypes
        show_freqNumLabel = true;
    end
    

    h_outer_top = (rowi - 1) * (h_textFreeze + h_colormap + h_deltay_colormap);
    h_inner_top = 0;
    h_inner_bottom = 0;
    h_outer_diff = h_colormap;
    if show_freezeType
        h_inner_top = h_inner_top + h_textFreeze;
        h_outer_diff = h_outer_diff + h_textFreeze;
    end
    if show_freqNumLabel
        h_inner_bottom = h_inner_bottom + h_textFreNum + h_textFreLabel;
        h_outer_diff = h_outer_diff + + h_textFreNum + h_textFreLabel;
    end
    h_outer_bottom = fig_height - (h_outer_top + h_outer_diff);
    
    w_outer_left = 0;
    w_outer_right = 0;
    w_inner_left = w_textpair;
    w_inner_right = w_textColorbar;
    
    % outer and inner margin
    subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
    subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
    
    titlename = freeType;
    titlename(1) = upper(titlename(1));
    titlename = strrep(titlename, 'Freeze', ' Freeze');
    plot_ciCohHistogram3(sigciCoh_flatten, chnPairNames, f_selected, titlename, ...
        'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ....
        'show_titlename', show_freezeType, 'show_xlabel', show_freqNumLabel, 'show_xticklabels', show_freqNumLabel, 'show_yticklabels', show_pairname);
    
    clear freeType;
end

%%% save
print(fullfile(savefolder, savefilename), '-dpng', '-r1000')
disp('saved figure')
close gcf

%%% compare cicoh event section
inputfolder_compare = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm3_fs500Hz_uNHP_histFTLagPhase');
compEvents = {'preMove', 'Early Reach'};


nrows = length(compEvents);
fig_height = h_colormap * nrows + h_deltay_colormap * (nrows -1) + h_textFreeze * nrows + (h_textFreNum + h_textFreLabel);
fig_width = w_textpair + w_colormap + w_textColorbar;
fig = figure('Position', [50 50 fig_width fig_height]);
set(fig, 'PaperUnits', 'points');
for rowi = 1 : nrows
    if rowi == 1
        ciCohPhasefile = ['Kitty ciCohPhasefile_8-40Hz_moderate_preMove_align2TargetOnset.mat'];
    else
        ciCohPhasefile = ['Kitty ciCohPhasefile_8-40Hz_moderate_earlyReach_align2ReachOnset.mat'];
    end
    
    load(fullfile(inputfolder_compare, ciCohPhasefile), 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
    
    % extract sigciCoh
    [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh);
    % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
    [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);
    
    
    
    show_freqNumLabel= false;
    show_pairname = true;
    show_freezeType = true;
    if rowi == nrows
        show_freqNumLabel = true;
    end
    

    h_outer_top = (rowi - 1) * (h_textFreeze + h_colormap + h_deltay_colormap);
    h_inner_top = 0;
    h_inner_bottom = 0;
    h_outer_diff = h_colormap;
    if show_freezeType
        h_inner_top = h_inner_top + h_textFreeze;
        h_outer_diff = h_outer_diff + h_textFreeze;
    end
    if show_freqNumLabel
        h_inner_bottom = h_inner_bottom + h_textFreNum + h_textFreLabel;
        h_outer_diff = h_outer_diff + + h_textFreNum + h_textFreLabel;
    end
    h_outer_bottom = fig_height - (h_outer_top + h_outer_diff);
    
    w_outer_left = 0;
    w_outer_right = 0;
    w_inner_left = w_textpair;
    w_inner_right = w_textColorbar;
    
    % outer and inner margin
    subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
    subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
    
    titlename = compEvents{rowi};
    plot_ciCohHistogram3(sigciCoh_flatten, chnPairNames, f_selected, titlename, ...
        'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ....
        'show_titlename', show_freezeType, 'show_xlabel', show_freqNumLabel, 'show_xticklabels', show_freqNumLabel, 'show_yticklabels', show_pairname);
end
print(fullfile(savefolder, [savefilename '-compared']), '-dpng', '-r1000')
disp('saved figure')
close gcf



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
    set(gca,'YTickLabel',chnPairNames,'fontsize',11, 'FontName', fontname, 'FontWeight', 'bold')
else 
    yticks([]);
end

if show_xlabel
    xlabel('Frequency (Hz)', 'fontsize', 12, 'FontName', fontname, 'FontWeight', 'bold')
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


