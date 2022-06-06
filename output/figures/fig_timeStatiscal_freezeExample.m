function fig_timeStatiscal_freezeExample()
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


[~, ~, ~, outputfolder] = exp_subfolders();

savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;

%% plot figure parameters in pixels
w_subplot = 360;
h_subplot_a = 300;
h_subplot_b1 = 200;
h_subplot_b2 = 200;

w_delta = 40;
h_delta_betweenab = 40;
h_delta_inb = 10;



w_ylabel = 60;
w_yticks = 20;
w_yright = 50;
h_xlabel = 30;
h_xticks = 20;
h_title = 40;

w_textab = 20;
h_textab = 20;


fig_width = w_textab +w_ylabel + w_yticks*2 + w_subplot * 2 + w_delta + w_yright;
fig_height = h_subplot_a + h_subplot_b1 + h_subplot_b2 + (h_xlabel + h_xticks) * 2 + h_delta_betweenab + h_delta_inb;

close all
fig = figure('Units', 'pixels', 'Position', [100 100 fig_width fig_height]);

%%%
h_top_a = 0;
h_top_b = h_subplot_a + (h_xticks) + h_delta_betweenab;
w_left_ab = w_textab;



%%% plot time statiscal
fig_timeStatiscal('fig', fig, 'h_top', h_top_a, 'w_left', w_left_ab, 'w_timePlot', w_subplot, 'h_timePlot', h_subplot_a, ...
                  'w_deltax_timePlots', w_delta, 'w_textTime', w_ylabel, 'w_textTimeNum', w_yticks, 'w_hold_yright', w_yright,...
                   'h_textAnimal', h_title, 'h_textCondLabel', h_xlabel);
              

              
fig_freezeExample('fig', fig, 'h_top', h_top_b, 'w_left', w_left_ab, 'w_subplotb', w_subplot, 'h_sublotb', h_subplot_b1, ...
                  'w_delta_b', w_delta  + h_xticks, 'h_delta_b', h_delta_inb, ...
                  'w_textSpeedFreq', w_ylabel, 'w_SpeedFreqTicks', w_yticks,  'w_colorbar', w_yright, 'holdplace_SpeedFreqTicks', true, ...
                  'h_textTime', h_xlabel, 'h_TimeTicks', h_xticks);



%%% plot text (a) and (b)
ta = annotation(fig, 'textbox', 'String', {'(a)'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
pos_left = 0;
pos_lower = fig_height - h_textab;
ta.Position = [pos_left pos_lower w_textab h_textab];

tb = annotation(fig, 'textbox', 'String', {'(b)'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
pos_lower = fig_height - (h_top_b + h_textab);
tb.Position = [pos_left pos_lower w_textab h_textab];

clear ta tb pos pos_left pos_lower


%%% save
print(fullfile(savefolder, savefilename), '-dpng', '-r1000')
print(gcf, fullfile(savefolder, savefilename), '-painters', '-depsc');
disp('saved figure')

function fig_timeStatiscal(varargin)
%   Inputs:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'h_top' - distance between the subplot top and fig top
%           'w_left' - distance between the subplot left and fig left
%           'w_timePlot' - width  for the time plot, default 360
%           'h_timePlot' - height for the time plot, default 300
%           'w_textTime' - width showing the 'Reachi Time (s)', default 40
%           'w_textTimeNum' - width showing the ylabel, default 40
%           'w_deltax_timePlots' - x distance between two time plots, default 50
%           'h_textAnimal' - height showing the animal, i.e. Animal J, default 40
%           'h_textCondLabel' - height showing the condition label, i.e. Normal, n = 176, default 30
%           'w_hold_yright' - width holding on the right, default 0 



% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'h_top', 0, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_left', 0, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_timePlot', 360, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_timePlot', 300, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_textTime', 40, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_textTimeNum', 10, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_deltax_timePlots', 50, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_textAnimal', 40, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_textCondLabel', 30, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_hold_yright', 0, @(x) assert(isscalar(x)&&isnumeric(x)));


parse(p,varargin{:});
fig = p.Results.fig;
h_top = p.Results.h_top;
w_left = p.Results.w_left;
w_timePlot = p.Results.w_timePlot;
h_timePlot = p.Results.h_timePlot;
w_textTime = p.Results.w_textTime;
w_textTimeNum = p.Results.w_textTimeNum;
w_deltax_timePlots = p.Results.w_deltax_timePlots;
h_textAnimal = p.Results.h_textAnimal;
h_textCondLabel = p.Results.h_textCondLabel;
w_hold_yright = p.Results.w_hold_yright;



%% Input 
[~, ~, pipelinefolder, ~] = exp_subfolders();
input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');




%% Code Start here
coli_reachonset = 2;
coli_touch = 3;

animals = {'Jo';'Kitty'};


conds_J = cond_cell_extract('Jo');
conds_K = cond_cell_extract('Kitty');

if isempty(fig)
    fig_width = w_textTime + w_textTimeNum*2 + w_timePlot * 2 + w_deltax_timePlots;
    fig_height = h_textAnimal + h_timePlot + h_textCondLabel;
    fig = figure('Position', [150 150 fig_width fig_height]);
    clear fig_width fig_height
end

pos = get(gcf, 'Position');
fig_width = pos(3);
fig_height = pos(4);
clear pos


%%% plot Reach time 
for ai = 1 : length(animals)
    animal = animals{ai};
    if ai == 1
        cond_cell = conds_J;
        input_folder = input_folder_J;
    else
        cond_cell = conds_K;
        input_folder = input_folder_K;
    end
    nconds = length(cond_cell);
    
    % extract t_reaction
    for ci = 1 : nconds
        pdcond  = cond_cell{ci};
        files = dir(fullfile(input_folder, ['*_' pdcond '_*.mat']));
        t_event = [];
        for fi = 1: length(files)
            load(fullfile(input_folder, files(fi).name), 'T_idxevent_lfp', 'fs_lfp');
            if ai == 1
                load(fullfile(input_folder, files(fi).name), 'goodTrials');
            else
                load(fullfile(input_folder, files(fi).name), 'selectedTrials');
                goodTrials = selectedTrials;
                clear selectedTrials
            end
            t_event = cat(1, t_event, T_idxevent_lfp{goodTrials, :}/ fs_lfp);
            clear T_idxevent_lfp fs_lfp
        end
        t_reaction.(pdcond) = t_event(:, coli_touch) - t_event(:, coli_reachonset);
        clear pdcond files t_event fi
    end
    
    show_timeLabel = false;
    if ai == 1
        show_timeLabel = true;
    end
    show_timeNum = true;
    show_animalLabel = true;
    show_condLabel = true;
    
    
    % w_out_left w_inner_left w_outer_right w_inner_right
    w_outer_left = w_left + w_textTime + (ai -1) * (w_textTimeNum + w_timePlot + w_deltax_timePlots);
    w_inner_left = 0;
    w_outer_diff = w_timePlot + w_hold_yright;
    w_inner_right = w_hold_yright;
    if show_timeLabel
        w_outer_left = w_outer_left - w_textTime;
        w_inner_left = w_inner_left + w_textTime;
        w_outer_diff = w_outer_diff + w_textTime;
    end
    if show_timeNum
        w_outer_diff = w_outer_diff + w_textTimeNum;
        w_inner_left = w_inner_left + w_textTimeNum;
    end
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);


    
    % h_outer_top h_inner_top h_outer_bottom h_inner_bottom
    h_outer_top = h_top;
    h_inner_top = 0;
    h_inner_bottom = 0;
    h_outer_diff = h_timePlot;
    if show_animalLabel
        h_inner_top = h_inner_top + h_textAnimal;
    end
    if show_condLabel
        h_inner_bottom = h_inner_bottom + h_textCondLabel;
        h_outer_diff = h_outer_diff + h_textCondLabel;
    end
    h_outer_bottom = fig_height -  (h_outer_top + h_outer_diff);
    
    % outer and inner margin
    subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
    subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
    
    plot_1timeStatiscal(t_reaction, 'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ...
        'show_xticklabels', show_condLabel, 'show_ylabel', show_timeLabel, 'show_yticklabels', show_timeNum, ...
        'titlename', ['Animal ' upper(animal(1))], 'FontName', 'Times New Roman')
    
    
    clear animal cond_cell nconds ci input_folder
    clear show_timeLabel show_animalLabel show_timeNum
    clear w_outer_left  w_inner_left w_outer_diff w_outer_right
    clear h_outer_top h_inner_top h_outer_bottom h_inner_bottom
    clear subplot_outerMargin subplot_innerposMargin
    clear t_reaction
end




function plot_1timeStatiscal(t_reaction, varargin)
%   Inputs:
%       t_reaction:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default [5 5 5 5]
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default [5 5 5 5]
%           'titlename' - title name for showing
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_ylabel' - show (true, default) or not show (false) ylabel


% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'innerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'titlename', '', @(x) assert(ischar(x)));
addParameter(p, 'fontname', 'Times New Roman', @ischar);
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
fig = p.Results.fig;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin; 
titlename = p.Results.titlename;
fontname = p.Results.fontname;
show_xticklabels = p.Results.show_xticklabels;
show_yticklabels = p.Results.show_yticklabels;
show_ylabel = p.Results.show_ylabel;



cond_cell = fieldnames(t_reaction);
nconds = length(cond_cell);

% grouping variables and ts for boxplot
ts = [];
gs =[];
for ci = 1 : length(cond_cell)
    pdcond  = cond_cell{ci};
    
    n = length(t_reaction.(pdcond));
    pdconds = repmat({[pdcond]}, n,1);
    gs = [gs; pdconds];
    
    ts = [ts; t_reaction.(pdcond)];
    
    clear pdcond n pdconds
end


%%% plot
if isempty(fig)
    fig = figure;
end

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
ax = axes(fig, 'Units', 'pixels');
set(gca, 'OuterPosition', outerpos, 'Position', innerpos)


% actual box plot
boxplot(ax, ts, gs); hold on
if ~isempty(titlename)
    title(titlename, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end
if show_ylabel
    ylabel('Reach Time (s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_xticklabels
    ax = gca;
    xtlabels = ax.XTickLabel;
    ax.XTickLabel = '';
    for xti = 1 : length(xtlabels)
        pdcond = xtlabels{xti};
        text2 = ['n=' num2str(length(t_reaction.(pdcond)))];
        newlabel = {pdcond; text2};
        text(xti, ax.YLim(1), sprintf('%s\n%s\n%s', newlabel{:,:}), 'horizontalalignment', 'center', ...
            'verticalalignment', ...
            'top', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', fontname);
    end
end



% significant part: wilcoxon rank sum test and plot
alpha = 0.5;
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
hs = [0.92 0.97 1.02];
for ci = 1 : nconds-1
    pdcondi  = cond_cell{ci};
    ti = t_reaction.(pdcondi);
    for cj = ci+1 : nconds
        pdcondj  = cond_cell{cj};
        tj = t_reaction.(pdcondj);
        p = ranksum(ti, tj);
        
        if p < alpha % plot sig * 
            if ci == 1 && cj == 2
                h = hs(1);
            elseif ci == 2 && cj == 3
                h = hs(2);
            elseif ci == 1 && cj == 3
                h = hs(3);
            end
            plot(xt([ci cj]), [1 1]*max(yt)*h, '-k',  mean(xt([ci cj])), max(yt) * (h+0.02), '*k')
            clear h
        end
        
        clear pdcondj tj 
    end
    clear cj
    clear pdcondi ti
end
clear ci

