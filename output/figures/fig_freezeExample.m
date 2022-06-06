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
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default [5 5 5 5]
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default [5 5 5 5]
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
addParameter(p, 'innerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
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

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
ax = axes(fig, 'Units', 'pixels');
set(gca, 'OuterPosition', outerpos, 'Position', innerpos)


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
    
    xtks = xticks(ax);
    xtklabs = xticklabels(ax);
    for tei = 1 : length(tevents_ma)
        tevent = tevents_ma(tei);
        xtks = [xtks tevent];
        xtklabs = [xtklabs; eveNames{tei}];
        clear tevent
    end
    [xtks, idxs]= sort(xtks);
    xtklabs = xtklabs(idxs);
    set(ax,'XTick', xtks, 'XTickLabel', xtklabs, 'XTickLabelRotation', 45, ...
        'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold');
    clear xtks xtklabs tevents_ma tei
end


% plot start and end freeze lines
ys = ylim;
ys = [ys(1) ys(2)/2];
for tfi = 1 : size(tFreezePhases, 1)
    plot([tFreezePhases(tfi, 1) tFreezePhases(tfi, 1)], ys, [colors4freezeline{tfi} '-'], 'LineWidth',1);
    plot([tFreezePhases(tfi, 2) tFreezePhases(tfi, 2)], ys, [colors4freezeline{tfi} '-'], 'LineWidth',1);
end
ar = annotation(fig,'doublearrow','Units','pixels');



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

