function fig2_timeStatiscal_freezeExample()
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
[~, ~, ~, outputfolder] = exp_subfolders();
[~, funcname, ~]= fileparts(codefilepath);


savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end


savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

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
                   'h_textAnimal', h_title, 'h_textCondLabel', h_xlabel, 'savefolder', savefolder);
              

              
fig_freezeExample('fig', fig, 'savefolder', savefolder, ...
                  'h_top', h_top_b, 'w_left', w_left_ab, 'w_subplotb', w_subplot, 'h_sublotb', h_subplot_b1, ...
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
%           'savefolder' - 



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
addParameter(p, 'savefolder', '', @isstr);


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
savefolder = p.Results.savefolder;



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
    
    
    ifig = figure('Position', [150 150 300 200]);
    plot_1timeStatiscal(t_reaction, 'fig', ifig, ...
        'show_xticklabels', show_condLabel, 'show_ylabel', show_timeLabel, 'show_yticklabels', show_timeNum, ...
        'titlename', ['Animal ' upper(animal(1))], 'FontName', 'Times New Roman')
    subfilename = ['time-' animal]; % 'Jo-mild-Bnormal-preMove'
    print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
    print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
    close(ifig)
    
    
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
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default []
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default []
%           'titlename' - title name for showing
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_ylabel' - show (true, default) or not show (false) ylabel


% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'innerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
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
ax = axes(fig, 'Units', 'pixels');
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
if ~isempty(outerposMargin) && ~~isempty(innerposMargin)
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
else
    pos = ax.Position;
    pos(2) = 30;
    pos(4) = fig_height - pos(2) - 20;
    set(gca, 'OuterPosition', [0 0 fig_width fig_height], 'Position', pos);
end



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


function fig_freezeExample(varargin)
%   Inputs:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'h_top' - distance between the subplot top and fig top, default 0
%           'w_left' - distance between the subplot left and fig left, default 0
%           'h_sublotb' - width  for the time plot, default 200
%           'w_subplotb' - height for the time plot, default 360
%           'h_delta_b' - width  for the time plot, default 20
%           'w_delta_b' - height for the time plot, default 10
%           'w_textSpeedFreq' - width  for the time plot, default 30
%           'w_SpeedFreqTicks' - height for the time plot, default 30
%           'w_colorbar' - width  for the time plot, default 40
%           'h_textTime' - height for the time plot, default 30
%           'h_TimeTicks' - width  for the time plot, default 20
%           'holdplace_SpeedFreqTicks' - hold place for SpeedFreqTicks (true) or not (false, default) 
%           'savefolder' - 



% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'h_top', 0, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_left', 0, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_sublotb', 360, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_subplotb', 300, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_delta_b', 20, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_delta_b', 10, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_textSpeedFreq', 30, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_SpeedFreqTicks', 30, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'w_colorbar', 40, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_textTime', 30, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'h_TimeTicks', 20, @(x) assert(isscalar(x)&&isnumeric(x)));
addParameter(p, 'holdplace_SpeedFreqTicks', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'savefolder', '.', @isstr);




parse(p,varargin{:});
fig = p.Results.fig;
h_top = p.Results.h_top;
w_left = p.Results.w_left;
h_sublotb = p.Results.h_sublotb;
w_subplotb = p.Results.w_subplotb;
h_delta_b = p.Results.h_delta_b;
w_delta_b = p.Results.w_delta_b;
w_textSpeedFreq = p.Results.w_textSpeedFreq;
w_SpeedFreqTicks = p.Results.w_SpeedFreqTicks;
w_colorbar = p.Results.w_colorbar;
h_textTime = p.Results.h_textTime;
h_TimeTicks = p.Results.h_TimeTicks;
holdplace_SpeedFreqTicks = p.Results.holdplace_SpeedFreqTicks;
savefolder = p.Results.savefolder;
 


%% Input 
[~, ~, pipelinefolder, ~] = exp_subfolders();

input_InitReachFreeze_file = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm3_fs500Hz_freezeSKTData_EpisodeExtract', 'Kitty_freezeEpisodes_moderate-tThesFreezeReach5s_20150408_bktdt2.mat');
input_ManiFreeze_file = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm3_fs500Hz_freezeSKTData_EpisodeExtract', 'Kitty_freezeEpisodes_moderate-tThesFreezeReach5s_20150408_bktdt2.mat');

tri_InitReachFreeze = 12;
tri_ManiFreeze = 4;

colors4freezeline_InitReachFreeze = {'k', 'b'};
colors4freezeline_ManiFreeze = {'r'};

freezeTypeInEpi_InitReachFreeze = {'freeze during init Move'; 'freeze during React-Reach'};
freezeTypeInEpi_ManiFreeze = {'freeze during Manipulation'};


%% Code Start here
ncols = 2;
nrows = 2;

if isempty(fig)
    fig_width = w_subplotb * 2 + w_delta_b;
    fig_height = h_sublotb * 2 + h_delta_b;
    fig = figure('Units', 'pixels', 'Position', [150 150 fig_width fig_height]);
    clear fig_width fig_height
end
pos = get(gcf, 'Position');
fig_width = pos(3);
fig_height = pos(4);
clear pos



for coli = 1: ncols
    show_TextSpeedFreq = false;
    show_Colorbar = false;
    show_SpeedFreqTicks = false;
    if coli == 1
        show_TextSpeedFreq = true;
        show_SpeedFreqTicks = true;
    end
    if coli == ncols
        show_Colorbar = true;
    end
    
    % w_out_left w_inner_left w_outer_right w_inner_right
    w_outer_left = w_left + w_textSpeedFreq + w_SpeedFreqTicks + (coli -1)* (w_subplotb + w_delta_b);
    w_inner_left = 0;
    w_inner_right = 0;
    w_outer_diff = w_subplotb;
    if show_TextSpeedFreq
        w_outer_left = w_outer_left - w_textSpeedFreq;
        w_inner_left = w_inner_left + w_textSpeedFreq;
        w_outer_diff = w_outer_diff + w_textSpeedFreq;
    end
    if show_SpeedFreqTicks || holdplace_SpeedFreqTicks
        w_outer_left = w_outer_left - w_SpeedFreqTicks;
        w_inner_left = w_inner_left + w_SpeedFreqTicks;
        w_outer_diff = w_outer_diff + w_SpeedFreqTicks;
    end   
    if show_Colorbar
        w_inner_right = w_inner_right + w_colorbar;
        w_outer_diff = w_outer_diff + w_colorbar;
    end
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);
    
    % load data
    if coli == 1
        input_Freeze_file = input_InitReachFreeze_file;
        tri = tri_InitReachFreeze;
        freezeTypeInEpis = freezeTypeInEpi_InitReachFreeze;
        colors4freezeline = colors4freezeline_InitReachFreeze;
    else
        input_Freeze_file = input_ManiFreeze_file;
        tri = tri_ManiFreeze;
        freezeTypeInEpis = freezeTypeInEpi_ManiFreeze;
        colors4freezeline = colors4freezeline_ManiFreeze;
    end
    load(input_Freeze_file, 'freezStruct', 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma');
    ma = smoothWspeed_trial{tri};
    ts = (1: length(ma))/fs_ma;
    tevents_ma = T_idxevent_ma{tri, :} / fs_ma;
    
    
    % find fei for each freezeTypeInEpi in freezeTypeInEpis
    freezEpisodes = freezStruct.freezEpisodes;  
    frzis = [];
    tFreezePhases = [];
    for fTi = 1: length(freezeTypeInEpis)
        freezeTypeInEpi = freezeTypeInEpis{fTi};
        
        % find freezEpisodes index fei
        frzi = 0;
        for fi = 1 : length(freezEpisodes)
            if(freezEpisodes{fi}.triali == tri && strcmp(freezEpisodes{fi}.freezeType, freezeTypeInEpi))
                frzi = fi;
                break;
            end
        end
        if frzi == 0
            disp(['Can not find freeze episode index for tri = ' num2str(tri) ':' freezeTyp]);
            return;
        end
        frzis = [frzis; frzi];
        
        tFreezePhases = [tFreezePhases; freezEpisodes{frzi}.freezeTPhaseS];
        
        clear freezeTypeInEpi fei
    end

    for rowi = 1 : nrows
        show_TextTime = false;
        show_TimeTicks = false;
        if rowi == nrows
            show_TextTime = true;
            show_TimeTicks = true;
        end
        
        % h_outer_top h_inner_top h_outer_bottom h_inner_bottom
        h_outer_top = h_top;
        h_inner_top = 0;
        h_inner_bottom = 0;
        h_outer_diff = h_sublotb;
        if show_TimeTicks
            h_inner_bottom = h_inner_bottom + h_TimeTicks;
            h_outer_diff = h_outer_diff + h_TimeTicks;
        end 
        if show_TextTime
            h_inner_bottom = h_inner_bottom + h_textTime;
            h_outer_diff = h_outer_diff + h_textTime;
        end
        h_outer_bottom = fig_height -  (h_outer_top + h_outer_diff);
        
        
        % outer and inner margin
        subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
        subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
        
        
        if rowi == 1
            plot_1freezeTrial(ma, ts, tevents_ma, tFreezePhases, colors4freezeline,...
                'fig', fig,...
                'show_xlabel', show_TextTime, 'show_xticklabels', show_TimeTicks, 'show_ylabel', show_TextSpeedFreq, 'show_yticklabels', show_SpeedFreqTicks, ...
                'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin)
            
            
            ifig = figure('Position', [150 150 400 300]);
            plot_1freezeTrial(ma, ts, tevents_ma, tFreezePhases, colors4freezeline,...
                'fig', ifig, 'show_xlabel', false)
            subfilename = ['freezTrial-' num2str(coli)]; 
            print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
            print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
            close(ifig)
        end
        
        
        clear show_TextTime show_TimeTicks
        clear h_outer_top h_inner_top h_outer_bottom h_inner_bottom
        clear subplot_outerMargin subplot_innerposMargin
        
    end
    
    clear show_TextSpeedFreq show_Colorbar show_SpeedFreqTicks
    clear w_outer_left w_outer_right w_inner_left w_inner_right w_outer_diff
    clear('freezStruct', 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma')
    clear ma ts tevents_ma freezEpisodes frzis tFreezePhases colors4freezeline
end


function plot_1freezeTrial(ma, ts, tevents_ma, tFreezePhases, colors4freezeline, varargin)
%   Inputs:
%       ma: 1d vector 
%       ts: 1d time points length same as ma
%       tevents_ma: event time points for ma data, 1 * nevents vector
%       tFreezePhase: start and end time point for all freeze Phases nfreezePhases * 2 [t_start t_end]
%       colors4freezeline:  color used for freeze line, length = nfreezePhases e.g. {'k', 'b', 'r'};
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default []
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default []
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_ylabel' - show (true, default) or not show (false) yticklabels
%           'show_xticklabels' - show (true, default) or not show (false) yticklabels
%           'show_eventLine' - show (true, default) or not show (false) event lines
%           'show_eventName' - show (true, default) or not show (false) event Names
%           'tlimit' - time show limit, default [] represents [min(ts) max(ts)]





% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'innerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'fontname', 'Times New Roman', @ischar);
addParameter(p, 'show_eventLine', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_eventName', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'tlimit', [], @(x) assert(isempty(x) || (isvector(x) && isnumeric(x) && length(x)==2)));




parse(p,varargin{:});
fig = p.Results.fig;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin;
show_xlabel = p.Results.show_xlabel;
show_xticklabels = p.Results.show_xticklabels;
show_ylabel = p.Results.show_ylabel;
show_yticklabels = p.Results.show_yticklabels;
fontname = p.Results.fontname;
show_eventLine = p.Results.show_eventLine;
show_eventName = p.Results.show_eventName;
tlimit = p.Results.tlimit;



if isempty(fig)
    fig = figure;
end

ax = axes(fig, 'Units', 'pixels');

% set axes position
if ~isempty(outerposMargin) && ~isempty(outerposMargin)
    fig_pos = fig.Position;
    fig_width = fig_pos(3);
    fig_height = fig_pos(4);
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
end



% plot MA
plot(ax, ts, ma, 'DisplayName','speed', 'LineWidth',1); 
hold on

% adjust xlim
if isempty(tlimit)
    tlimit = [min(ts) max(ts)];
end
xlim(tlimit);

% plot threshold
speedThres_Move = 30;
plot(xlim, [speedThres_Move speedThres_Move], 'b-.');
%text(pi,0,'\leftarrow sin(\pi)')


% plot event line
if show_eventLine
    for tei = 1 : length(tevents_ma)
        tevent = tevents_ma(tei);
        plot([tevent tevent], ylim, '--');
        clear tevent
    end
end


% plot event Name
if show_eventName
    eveNames = {'CueOnset', 'ReachOnset', 'Touch', 'ReturnOnset', 'Mouth'};
    
    xtks = [];
    xtklabs = [];
    for tei = 1 : length(tevents_ma)
        tevent = tevents_ma(tei);
        xtks = [xtks tevent];
        xtklabs = [xtklabs; {eveNames{tei}}];
        clear tevent
    end
    [xtks, idxs]= sort(xtks);
    xtklabs = xtklabs(idxs);
    set(ax,'XTick', xtks, 'XTickLabel', xtklabs, 'XTickLabelRotation', 45, ...
        'FontName', 'Times New Roman');
    clear xtks xtklabs tei
end


% plot start and end freeze lines
ys = ylim;
ys = [ys(1) ys(2)/2];
for tfi = 1 : size(tFreezePhases, 1)
    plot([tFreezePhases(tfi, 1) tFreezePhases(tfi, 1)], ys, [colors4freezeline{tfi} '-'], 'LineWidth',1);
    plot([tFreezePhases(tfi, 2) tFreezePhases(tfi, 2)], ys, [colors4freezeline{tfi} '-'], 'LineWidth',1);
end

% plot t_touch and t_returnonset as ManiFreeze
eveNames = {'CueOnset', 'ReachOnset', 'Touch', 'ReturnOnset', 'Mouth'};
t_touch = tevents_ma(cellfun(@(x) strcmp(x,'Touch'), eveNames));
t_returnonset = tevents_ma(cellfun(@(x) strcmp(x,'ReturnOnset'), eveNames));
plot([t_touch t_touch], ys, ['r-'], 'LineWidth',1);
plot([t_returnonset t_returnonset], ys, ['r-'], 'LineWidth',1);


if show_xticklabels
else
     xticks([]);
end

if show_xlabel
    xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_yticklabels
    set(ax.YAxis,'FontSize', 11, 'FontWeight', 'normal', 'FontName', fontname);
else
     yticks([]);
end

if show_ylabel
    ylabel('Speed', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end



