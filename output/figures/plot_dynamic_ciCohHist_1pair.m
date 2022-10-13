function plot_dynamic_ciCohHist_1pair(ciCoh, f_selected, t_selected, time0name,...
    varargin)
%
%   Inputs:

% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'cbarTicks', [], @(x) assert(isvector(x) && isnumeric(x)));
addParameter(p, 'cbarStr', '', @isstr);
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_titlename', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'titlename', '', @isstr);

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig = p.Results.fig;
cbarTicks = p.Results.cbarTicks;
cbarStr = p.Results.cbarStr;
show_ylabel = p.Results.show_ylabel;
show_xlabel = p.Results.show_xlabel;
show_titlename = p.Results.show_titlename;
show_colorbar = p.Results.show_colorbar;
titlename = p.Results.titlename;


% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



%%%  plot 

if isempty(fig)
    width = p.Results.width;
    height = p.Results.height;
    fig = figure('Position', [50 50 width height]);
    set(fig, 'PaperUnits', 'points');
    clear width height
end

ax = axes(fig, 'Units', 'pixels');

imagesc(ax, t_selected, f_selected, ciCoh); hold on;
colormap(jet)
set(gca,'YDir','normal')
if ~isempty(histClim)
    set(gca,'CLim', histClim)
end


% set time0name
xtklabels = xticklabels();
xtklabels{cellfun(@(x) strcmp(x, '0'), xtklabels)} = time0name;
xticklabels(xtklabels);

% time0name line
plot([0 0], ylim, '--r', 'LineWidth',2)


c = colorbar;

if ~isempty(cbarStr)
    c.Label.String = cbarStr;
end
if ~isempty(cbarTicks)
    set(c, 'Ticks', cbarTicks, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roma')
end
c.Visible = 'off';


%%% show inf
if show_ylabel
    ylabel('Frequence (Hz)','FontName','Times New Roma')
end

if show_xlabel
    xlabel('time/s','FontName','Times New Roma')
end

if show_titlename && ~isempty(titlename)
    title(titlename, 'FontSize', 10, 'FontWeight', 'normal')
end

if show_colorbar
    c.Visible = 'on';
end