function plot_ciCohHistogram(ciCoh_flatten, chnPairNames, f_selected, titlename, varargin)
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
%           'fig_left' - figure position left (default 50)
%           'fig_bottom' - figure position bottom (default 50)
%           'fig_width' - figure position width (default 1200)
%           'fig_height' - figure position height (default 60)
%           'cbarTicks' - vector, default [0 0.5 1]
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_titlename' - show (true, default) or not show (false) titlename
%           'show_colorbar' - show (true, default) or not show (false) colorbar
%



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [0 1], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'fig_left', 50, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_bottom', 50, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_width', 1000, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_height', 250, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'cbarTicks', [0 0.5 1], @(x) assert(isvector(x) && isnumeric(x)));
addParameter(p, 'cbarStr', 'ciCoh', @isstr);
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_titlename', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'margin_inner', [80 30 15 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==4));
addParameter(p, 'margin_outer', [5 5 5 5], @(x) assert(isnumeric(x) && isvector(x) && length(x)==4));

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig_left = p.Results.fig_left;
fig_bottom = p.Results.fig_bottom;
fig_width = p.Results.fig_width;
fig_height = p.Results.fig_height;
cbarTicks = p.Results.cbarTicks;
cbarStr = p.Results.cbarStr;
show_xticklabels = p.Results.show_xticklabels;
show_yticklabels = p.Results.show_yticklabels;
show_xlabel = p.Results.show_xlabel;
show_titlename = p.Results.show_titlename;
show_colorbar = p.Results.show_colorbar;
margin_inner = p.Results.margin_inner; % left, top, right, bottom of distance between outer and inner position
margin_outer = p.Results.margin_outer; % left, top, right, bottom margin of outer position


% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% plot
figure;
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
imagesc(ciCoh_flatten)
colormap(jet)


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
    chnPairNames{ci} = chnPair;
    chnPair_prev = chnPair;
    
    clear s_stn s_gp chnPair
end


% set axes position
outpos = [margin_outer(1) margin_outer(4) fig_width-margin_outer(1)-margin_outer(3) fig_height-margin_outer(2)-margin_outer(4)];
innerpos = [outpos(1)+margin_inner(1) outpos(2)+margin_inner(4) outpos(3)-margin_inner(1)-margin_inner(3) outpos(4)-margin_inner(2)-margin_inner(4)];
set(gca, 'Units', 'pixels');
set(gca, 'OuterPosition', outpos, 'Position', innerpos)
set(gca,'CLim', histClim)
c = colorbar;
c.Label.String = cbarStr;
if ~isempty(cbarTicks)
    c.Ticks = cbarTicks;
end
c.Visible = 'off';

%%% show inf
[npairs, nf] = size(ciCoh_flatten);
if show_xticklabels
    xticks([1:nf])
    xticklabels(round(f_selected))
else
    xticks([]);
    margin_inner(4) = margin_inner(4)-10;
end

if show_yticklabels
    yticks([1:npairs]);
    set(gca,'YTickLabel',chnPairNames,'fontsize',10,'FontWeight','normal')
else 
    yticks([]);
    margin_inner(1) = 0;
end

if show_xlabel
    xlabel('Frequence (Hz)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName','Times New Roma')
else
    margin_inner(4) = margin_inner(4)-30;
end
if show_titlename
    title(titlename, 'FontSize', 10, 'FontWeight', 'normal')
else
    margin_inner(2) = 0;
end

if show_colorbar
c.Visible = 'on';
end


