function fig5_freezeExample()
%   
%   Usage:
%       fig5_freezeExample()
%



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
fig_freezeExample('pos_ifig', [150 150 500 250],...
    'savefolder', savefolder, 'copy2folder', aisavefolder);



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

tri_Mani = 4;
tri_InitReach = 12;
colors4freezeline = {'k', 'b', 'c'};

% load data
load(input_Freeze_file, 'freezStruct', 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma');


% find tFreezePhases for each freezeTypeInEpi in freezeTypeInEpis
tri = tri_InitReach;
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

tFreezePhases_InitReach = tFreezePhases;
clear tFreezePhases


tri = tri_Mani;
freezEpisodes = freezStruct.freezEpisodes;
tFreezePhases = [];
freezeTypeInEpis = {'freeze during Manipulation'};
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

tFreezePhases_Mani = tFreezePhases;
clear tFreezePhases


% combined freeze phases from two trials
idx_Touch_inInitReach = T_idxevent_ma{tri_InitReach, 3};
ma_befTouch = smoothWspeed_trial{tri_InitReach}(1:idx_Touch_inInitReach);
idx_Touch_inMani = T_idxevent_ma{tri_Mani, 3};
ma_aftTouch = smoothWspeed_trial{tri_Mani}(idx_Touch_inMani+1:end);
ma = [ma_befTouch; ma_aftTouch];

idxevent_ma_befTouch = T_idxevent_ma{tri_InitReach, 1:3};
idxevent_ma_aftTouch = T_idxevent_ma{tri_Mani, 4:5} - idx_Touch_inMani + idx_Touch_inInitReach;
tevents_ma = [idxevent_ma_befTouch idxevent_ma_aftTouch] / fs_ma;
ts = (1: length(ma))/fs_ma;
tFreezePhases = [tFreezePhases_InitReach; tFreezePhases_Mani + ( -idx_Touch_inMani + idx_Touch_inInitReach)/fs_ma];
clear idx_Touch_inInitReach ma_befTouch idx_Touch_inMani ma_aftTouch
clear idxevent_ma_befTouch idxevent_ma_aftTouch

% + 2s for freeze onset of Init and Mani Freeze 
tFreezePhases(1, 1) = tFreezePhases(1, 1) +2;
tFreezePhases(3, 1) = tFreezePhases(3, 1) +2;

% plot
ifig = figure('Position', pos_ifig);
plot_1freezeTrial(ma, ts, tevents_ma, tFreezePhases, colors4freezeline,...
    'fig', ifig, 'show_xlabel', false, 'plotFreezeLines', true)
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
%       tFreezePhase: start and end time point for all freeze Phases nfreezePhases * 3 [t_start t_end]
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
        plot([tevent tevent], [0 speedThres_Move], '--');
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
        plot([tFreezePhases(tfi, 1) tFreezePhases(tfi, 1)], ys, 'k-', 'LineWidth',1);
        plot([tFreezePhases(tfi, 2) tFreezePhases(tfi, 2)], ys, 'k-', 'LineWidth',1);
    end

    % add cueonset+2 and touch + 2 xticklabels
    xtks = get(ax, 'XTick');
    xtklabs = get(ax, 'XTickLabel');
    xtks = [xtks tFreezePhases(1, 1)];
    xtklabs = [xtklabs; 'CueOnset+2'];
    xtks = [xtks tFreezePhases(3, 1)];
    xtklabs = [xtklabs; 'Touch+2'];
    [xtks, idxs]= sort(xtks);
    xtklabs = xtklabs(idxs);
    set(ax,'XTick', xtks, 'XTickLabel', xtklabs, 'XTickLabelRotation', 45, ...
        'FontName', 'Times New Roman');

    clear xtks xtklabs
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
