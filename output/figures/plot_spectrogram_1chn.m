function plot_spectrogram_1chn(psd_plot, freqs_plot, times_plot, align2, clim, varargin)
%
%   Inputs:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default []
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default []
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_ylabel' - show (true, default) or not show (false) ylabel
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_colorbar' - show (true, default) or not show (false) colorbar
%



% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'outerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'innerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
fig = p.Results.fig;
outerposMargin = p.Results.outerposMargin; 
innerposMargin = p.Results.innerposMargin;
show_xlabel = p.Results.show_xlabel;
show_xticklabels = p.Results.show_xticklabels;
show_ylabel = p.Results.show_ylabel;
show_yticklabels = p.Results.show_yticklabels;
show_colorbar = p.Results.show_colorbar;



if isempty(fig)
    fig = figure;
end

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
ax = axes(fig, 'Units', 'pixels');
if ~isempty(outerposMargin) && ~isempty(innerposMargin)
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
end


% plot
imagesc(ax, times_plot, freqs_plot, psd_plot)
set(ax,'YDir','normal', 'CLim', clim)
colormap(jet)
c = colorbar;
c.Visible = 'off';

%%% show inf
if show_xlabel
  xlabel('time (s)', 'FontSize', 12, 'FontWeight', 'bold')
end

if show_xticklabels
    xtls = xticklabels(ax);
    xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
    xticklabels(ax, xtls)
else
    xticks([]);
end

if show_ylabel
  ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
end

if show_yticklabels
else
    yticks([]);
end

if show_colorbar
  c.Visible = 'on';
end

end