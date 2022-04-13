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

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig_left = p.Results.fig_left;
fig_bottom = p.Results.fig_bottom;
fig_width = p.Results.fig_width;
fig_height = p.Results.fig_height;
cbarTicks = p.Results.cbarTicks;
cbarStr = p.Results.cbarStr;


% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% plot
dispix_outinpos = [80 30 15 40]; % left, top, right, bottom of distance between outer and inner position
margin_outpos = [5 5 5 5]; % left, top, right, bottom margin of outer position
outpos = [margin_outpos(1) margin_outpos(4) fig_width-margin_outpos(1)-margin_outpos(3) fig_height-margin_outpos(2)-margin_outpos(4)];
innerpos = [outpos(1)+dispix_outinpos(1) outpos(2)+dispix_outinpos(4) outpos(3)-dispix_outinpos(1)-dispix_outinpos(3) outpos(4)-dispix_outinpos(2)-dispix_outinpos(4)];
figure;
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
imagesc(ciCoh_flatten)
colormap(jet)
set(gca, 'Units', 'pixels');
set(gca, 'OuterPosition', outpos, 'Position', innerpos)
[npairs, nf] = size(ciCoh_flatten);
xticks([1:nf])
xticklabels(round(f_selected))
yticks([1:npairs]);
set(gca,'YTickLabel',chnPairNames,'fontsize',10,'FontWeight','normal')
xlabel('freqs/Hz')
title(titlename, 'FontSize', 10, 'FontWeight', 'normal')
set(gca,'CLim', histClim)
c = colorbar;
c.Label.String = cbarStr;
if ~isempty(cbarTicks)
    c.Ticks = cbarTicks;
end

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
    chnPair_prev = chnPair;
    
    clear s_stn s_gp chnPair
end