function fig1_reachTimeStatiscal(varargin)
%   
%   Usage:
%       fig1_reachTimeStatiscal('plot_timeStatiscal', false, 'plot_reachTimehist', true)
%       fig1_reachTimeStatiscal('plot_timeStatiscal', true, 'pos_ifig_reachTimeStatiscal', [150 150 300 250])
%
%   Inputs:
%
%       Name-Value:
%           'pos_ifig_reachTimeStatiscal' - position and size of the reachtime statistical figure [left bottom fig_width fig_height], default [150 150 300 250]
%           'pos_ifig_reachTimehist' - position and size of the reachtime hist figure [left bottom fig_width fig_height], default [150 150 300 250]
%           'plot_timeStatiscal' - tag plotting timeStatiscal (default true)
%           'plot_reachTimehist' - tag plotting reach time histogram (default false)


% parse params
p = inputParser;
addParameter(p, 'pos_ifig_reachTimeStatiscal', [150 150 300 250], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'pos_ifig_reachTimehist', [150 150 300 250], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'plot_timeStatiscal', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plot_reachTimehist', false, @(x) assert(islogical(x) && isscalar(x)));


parse(p,varargin{:});
pos_ifig_reachTimeStatiscal = p.Results.pos_ifig_reachTimeStatiscal;
pos_ifig_reachTimehist = p.Results.pos_ifig_reachTimehist;
plot_timeStatiscal = p.Results.plot_timeStatiscal;
plot_reachTimehist = p.Results.plot_reachTimehist;


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

aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', funcname);
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

%% Code start Here
if plot_timeStatiscal
    fig_timeStatiscal('pos_ifig', pos_ifig_reachTimeStatiscal, ...
        'savefolder', savefolder, 'copy2folder', aisavefolder);
end



if plot_reachTimehist
    fig_reachTimehist('pos_ifig', pos_ifig_reachTimehist,...
        'savefolder', savefolder, 'copy2folder', aisavefolder);
end


function fig_reachTimehist(varargin)
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the figure [left bottom fig_width fig_height], default [150 150 400 300]
%
%           'savefolder'
%           'copy2folder'

% parse params
p = inputParser;
addParameter(p, 'savefolder', '.', @ischar);
addParameter(p, 'pos_ifig', [150 150 400 300], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'copy2folder', '', @ischar);

parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;
savefolder = p.Results.savefolder;
copy2folder = p.Results.copy2folder;



% Input
[~, ~, pipelinefolder, ~] = exp_subfolders();
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');

conds_K = cond_cell_extract('Kitty');

coli_reachonset = 2;
coli_touch = 3;

animals = {'Kitty'};


%  Reach time histogram
for ai = 1 : length(animals)
    animal = animals{ai};
    if strcmpi(animal, 'Kitty')
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
            if strcmpi(animal, 'Jo')
                load(fullfile(input_folder, files(fi).name), 'goodTrials');
            else
                load(fullfile(input_folder, files(fi).name), 'selectedTrials');
                goodTrials = selectedTrials;
                clear selectedTrials
            end
            t_event = cat(1, t_event, T_idxevent_lfp{goodTrials, :}/ fs_lfp);
            clear T_idxevent_lfp fs_lfp goodTrials
        end
        t_reach.(pdcond) = t_event(:, coli_touch) - t_event(:, coli_reachonset);
        clear pdcond files t_event fi
    end


    ifig = figure('Position', pos_ifig);
    ax = axes(ifig);
    histogram(ax, t_reach.moderate, 25);
    title(['reach time histogram, ntrials = ' num2str(length(t_reach.moderate))])

    xlimits = xlim();
    tmp1 = round(xlimits(1)*10)/10;
    tmp2 = round(xlimits(2)*10)/10;
    xticks([tmp1 : 0.2: tmp2 ])
    
    subfilename = ['reachTimeHist-' animal]; 
    savefile = fullfile(savefolder, subfilename);
    print(ifig, savefile, '-painters', '-depsc')
    print(ifig, savefile, '-dpng', '-r1000')
    if ~isempty(copy2folder)
        print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
    end
    close(ifig)


    clear animal cond_cell input_folder nconds ci 
    clear t_reach
    clear show_timeLabel show_timeNum show_condLabel
    clear ifig subfilename savefile
end



function fig_timeStatiscal(varargin)
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the figure [left bottom fig_width fig_height], default [150 150 400 300]
%
%           'savefolder'
%           'copy2folder'

% parse params
p = inputParser;
addParameter(p, 'savefolder', '.', @ischar);
addParameter(p, 'pos_ifig', [150 150 400 300], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'copy2folder', '', @ischar);

parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;
savefolder = p.Results.savefolder;
copy2folder = p.Results.copy2folder;



% Input
[~, ~, pipelinefolder, ~] = exp_subfolders();
input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');




% Code Start here
coli_reachonset = 2;
coli_touch = 3;

animals = {'Jo';'Kitty'};


conds_J = cond_cell_extract('Jo');
conds_K = cond_cell_extract('Kitty');



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
            clear T_idxevent_lfp fs_lfp goodTrials
        end
        t_reach.(pdcond) = t_event(:, coli_touch) - t_event(:, coli_reachonset);
        clear pdcond files t_event fi
    end

    show_timeLabel = false;
    if ai == 1
        show_timeLabel = true;
    end
    show_timeNum = true;
    show_condLabel = true;

    ifig = figure('Position', pos_ifig);
    plot_1timeStatiscal(t_reach, 'fig', ifig, ...
        'show_xticklabels', show_condLabel, 'show_ylabel', show_timeLabel, 'show_yticklabels', show_timeNum, ...
        'titlename', ['Animal ' upper(animal(1))], 'FontName', 'Times New Roman')
    subfilename = ['time-' animal]; % 'Jo-mild-Bnormal-preMove'
    savefile = fullfile(savefolder, subfilename);
    print(ifig, savefile, '-painters', '-depsc')
    print(ifig, savefile, '-dpng', '-r1000')
    if ~isempty(copy2folder)
        print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
    end
    close(ifig)


    clear animal cond_cell input_folder nconds ci 
    clear t_reach
    clear show_timeLabel show_timeNum show_condLabel
    clear ifig subfilename savefile
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
hs = [0.88 0.92 0.97];
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
%           'pos_ifig' - position and size of the figure [left bottom fig_width fig_height], default [150 150 400 300]
%
%           'savefolder'
%           'copy2folder'

% parse params
p = inputParser;
addParameter(p, 'savefolder', '.', @ischar);
addParameter(p, 'pos_ifig', [150 150 400 300], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'copy2folder', '', @ischar);

parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;
savefolder = p.Results.savefolder;
copy2folder= p.Results.copy2folder;

% Input 
[~, ~, pipelinefolder, ~] = exp_subfolders();
input_Freeze_file = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm3_fs500Hz_freezeSKTData_EpisodeExtract', 'Kitty_freezeEpisodes_moderate-tThesFreezeReach5s_20150408_bktdt2.mat');

tri = 12;
colors4freezeline = {'k', 'b'};

% load data
load(input_Freeze_file, 'freezStruct', 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma');

ma = smoothWspeed_trial{tri};
ts = (1: length(ma))/fs_ma;
tevents_ma = T_idxevent_ma{tri, :} / fs_ma;

% find tFreezePhases for each freezeTypeInEpi in freezeTypeInEpis
freezEpisodes = freezStruct.freezEpisodes;
tFreezePhases = [];
freezeTypeInEpis = {'freeze during init Move'; 'freeze during React-Reach'};
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

    tFreezePhases = [tFreezePhases; freezEpisodes{frzi}.freezeTPhaseS];

    clear freezeTypeInEpi fei
end


ifig = figure('Position', pos_ifig);
plot_1freezeTrial(ma, ts, tevents_ma, tFreezePhases, colors4freezeline,...
    'fig', ifig, 'show_xlabel', false)
subfilename = 'freezTrial';
print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
if ~isempty(copy2folder)
    print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
end
close(ifig)

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
%           'plotFreezeLines' - plot freeze lines or not (defalut false)




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
addParameter(p, 'plotFreezeLines', false, @(x) assert(islogical(x) && isscalar(x)));




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
plotFreezeLines = p.Results.plotFreezeLines;



if isempty(fig)
    fig = figure;
end

ax = axes(fig, 'Units', 'pixels');

% set axes position
if ~isempty(outerposMargin) && ~isempty(innerposMargin)
    fig_pos = fig.Position;
    fig_width = fig_pos(3);
    fig_height = fig_pos(4);
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
end



% plot MA & threshold
plot(ax, ts, ma, 'DisplayName','speed', 'LineWidth',1); 
hold on

% adjust xlim
if isempty(tlimit)
    tlimit = [min(ts) max(ts)];
end
xlim(tlimit);

% plot threshold
speedThres_Move = 30;
plot(xlim, [speedThres_Move speedThres_Move], 'r-.');
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

if plotFreezeLines

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
end




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
    ylabel('Speed of Wrist (cm/s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end
