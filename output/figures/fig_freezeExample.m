function fig_freezeExample()


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
input_Freeze_file = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm3_fs500Hz_freezeSKTData_EpisodeExtract', 'Kitty_freezeEpisodes_moderate-tThesFreezeReach5s_20150408_bktdt2.mat');



savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;


%% Code Start here
load(input_Freeze_file, 'freezStruct', 'fs_ma', 'smoothWspeed_trial', 'T_idxevent_ma');



freezEpisodes = freezStruct.freezEpisodes;

%%% find freezEpisodes index fei
tri = 4;
freezeTypeInEpi = 'freeze during Manipulation';
fei = 0;
for fi = 1 : length(freezEpisodes)
    if(freezEpisodes{fi}.triali == tri && strcmp(freezEpisodes{fi}.freezeType, freezeTypeInEpi))
        fei = fi;
        break;
    end
end
if fei == 0
    disp(['Can not find freeze episode index for tri = ' num2str(tri) ':' freezeTyp]);
    return;
end
fei_mani = fei;
tri_mani = tri;
clear tri freezeTypeInEpi

%%% plot
h_sublotb = 200;
w_subplotb = 360;


h_delta_b = 20;
w_delta_b = 10;

w_textSpeedFreq = 30;
w_SpeedFreqTicks = 30;
w_colorbar = 40;

h_textTime = 30;
h_TimeTicks = 20;

fig_width = w_subplotb * 2 + w_delta_b;
fig_height = h_sublotb * 2 + h_delta_b;
fig = figure('Units', 'pixels', 'Position', [150 150 fig_width fig_height]);

ncols = 2;
nrows = 2;

tris = [12, 4];


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
    w_outer_left = w_textSpeedFreq + (coli -1)* (w_subplotb + w_delta_b);
    w_inner_left = 0;
    w_inner_right = 0;
    w_outer_diff = w_subplotb;
    if show_TextSpeedFreq
        w_outer_left = w_outer_left - w_textSpeedFreq;
        w_inner_left = w_inner_left + w_textSpeedFreq;
        w_outer_diff = w_outer_diff + w_textSpeedFreq;
    end
    if show_SpeedFreqTicks
        w_outer_left = w_outer_left - w_SpeedFreqTicks;
        w_inner_left = w_inner_left + w_SpeedFreqTicks;
        w_outer_diff = w_outer_diff + w_SpeedFreqTicks;
    end
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);
    if show_Colorbar
        w_inner_right = w_inner_right + w_colorbar;
    end
    
    tri = tris(coli);
    ma = smoothWspeed_trial{tri};
    ts = (1: length(ma))/fs_ma;
    
    for rowi = 1 : nrows
        show_TextTime = false;
        show_TimeTicks = false;
        if rowi == nrows
            show_TextTime = true;
            show_TimeTicks = true;
        end
        
        % h_outer_top h_inner_top h_outer_bottom h_inner_bottom
        h_outer_top = 0;
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
            plot_1freezeTrial(ma, ts, 'fig', fig,...
            'show_xlabel', show_TextTime, 'show_xticklabels', show_TimeTicks, 'show_ylabel', show_TextSpeedFreq, 'show_yticklabels', show_SpeedFreqTicks, ...
            'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin)
        end
        
        
        clear show_TextTime show_TimeTicks
        clear h_outer_top h_inner_top h_outer_bottom h_inner_bottom
        clear subplot_outerMargin subplot_innerposMargin
        
    end
    
    clear show_TextSpeedFreq show_Colorbar
    clear w_outer_left w_outer_right w_inner_left w_inner_right w_outer_diff
end

freezeType = 'Manipulate Freeze';








function plot_1freezeTrial(ma, ts, varargin)
%   Inputs:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default [5 5 5 5]
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default [5 5 5 5]
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_ylabel' - show (true, default) or not show (false) yticklabels
%           'show_xticklabels' - show (true, default) or not show (false) yticklabels


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



parse(p,varargin{:});
fig = p.Results.fig;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin;
show_xlabel = p.Results.show_xlabel;
show_xticklabels = p.Results.show_xticklabels;
show_ylabel = p.Results.show_ylabel;
show_yticklabels = p.Results.show_yticklabels;
fontname = p.Results.fontname;



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
hp = plot(ax, ts, ma, 'DisplayName','speed'); 
hold on

if show_xlabel
    xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_xticklabels
else
     xticks([]);
end

if show_ylabel
    ylabel('Speed', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_yticklabels
else
     yticks([]);
end

