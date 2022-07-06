function fig1_normalTrial_SpeedExample(varargin)
%   
%   Usage:
%       fig1_normalTrial_SpeedExample('pos_ifig', [150 150 300 250])
%
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the reachtime statistical figure [left bottom fig_width fig_height], default [150 150 300 250]
%  


% parse params
p = inputParser;
addParameter(p, 'pos_ifig', [150 150 300 250], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));


parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;


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
fig_normaltrialSpeed_Example('pos_ifig', pos_ifig,...
    'savefolder', savefolder, 'copy2folder', aisavefolder);



function fig_normaltrialSpeed_Example(varargin)
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
input_file = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'm2_segSKTData_SelectTrials_goodReach', 'Kitty_TrialsWMarkers_moderate_20150408_bktdt2.mat');

tri = 2;


% load data
load(input_file, 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma');

ma = smoothWspeed_trial{tri};
ts = [1 : length(ma)] / fs_ma;
tevents_ma = T_idxevent_ma{tri, :} / fs_ma;


% plot
ifig = figure('Position', pos_ifig);
plot_1freezeTrial(ma, ts, tevents_ma, ...
    'fig', ifig, 'show_xlabel', false)
subfilename = 'normalTrial_Speed_Example';
print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
if ~isempty(copy2folder)
    print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
end
close(ifig)

function plot_1freezeTrial(ma, ts, tevents_ma, varargin)
%   Inputs:
%       ma: 1d vector 
%       ts: 1d time points length same as ma
%       tevents_ma: event time points for ma data, 1 * nevents vector
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
    set(ax,'XTick', xtks, 'XTickLabel', xtklabs, 'XTickLabelRotation', 15, ...
        'FontName', 'Times New Roman');
    clear xtks xtklabs tei
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
