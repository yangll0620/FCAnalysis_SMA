function fig_timeStatiscal()


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
input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');


savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;


%% plot figure parameters in pixals
w_timePlot = 360; % width  for the time plot
h_timePlot = 300; % height for the time plot

w_textTime = 40; % width showing the 'Reachi Time (s)'
w_textTimeNum = 10; % width showing the ylabel
w_deltax_timePlots = 50; % x distance between two time plots


h_deltay = 15; % y distance between two plots
h_textAnimal = 40; % height showing the animal, i.e. Animal J
h_textCondLabel = 30; % height showing the condition label, i.e. Normal, n = 176

fontSize_textAnimal = 12;
fontname = 'Times New Roman';

%% Code Start here
coli_reachonset = 2;
coli_touch = 3;

animals = {'Jo';'Kitty'};


conds_J = cond_cell_extract('Jo');
conds_K = cond_cell_extract('Kitty');


fig_width = w_textTime + w_textTimeNum*2 + w_timePlot * 2 + w_deltax_timePlots;
fig_height = h_textAnimal + h_timePlot + h_textCondLabel;

fig = figure('Units', 'pixels', 'Position', [100 100 fig_width fig_height]);

%%% plot text (a) and (b)
w_text = 40;
h_text = 20;
t1 = annotation(fig, 'textbox', 'String', {'(a)'}, 'LineStyle', 'none', 'Units', 'pixels', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname);
pos_left = 0;
pos_lower = fig_height - h_text;
t1.Position = [pos_left pos_lower w_text h_text];
clear t1 pos pos_left pos_lower

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
    w_outer_left = w_textTime + (ai -1) * (w_textTimeNum + w_timePlot + w_deltax_timePlots);
    w_inner_left = 0;
    w_outer_diff = w_timePlot;
    w_inner_right = 0;
    if show_timeLabel
        w_outer_left = w_outer_left - w_textTime;
        w_inner_left = w_inner_left + w_textTime;
        w_outer_diff = w_outer_diff + w_textTime;
    end
    if show_timeNum
        w_outer_diff = w_outer_diff + w_textTimeNum;
    end
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);
    if ai == 2
        w_outer_right = w_outer_right + 10;
    end

    
    % h_outer_top h_inner_top h_outer_bottom h_inner_bottom
    h_outer_top = 0;
    h_inner_top = 0;
    h_outer_diff = h_timePlot;
    if show_animalLabel
        h_inner_top = h_inner_top + h_textAnimal;
    end
    if show_condLabel
        h_inner_bottom = h_textCondLabel;
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
%       cond_cell:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default [5 5 5 5]
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default [5 5 5 5]
%           'titlename' - title name for showing


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
    ylabel('Reach Time/s', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
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

