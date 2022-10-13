function plot_ciCohHist_dynamic(ciCoh, f_selected, t_selected, titlename, time0name, histClim, varargin)
%
%   Inputs:



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% plot
figure;
imagesc(t_selected, f_selected, ciCoh); hold on;
colormap(jet)
set(gca,'YDir','normal')
set(gca,'CLim', histClim)

ylabel('Frequence (Hz)','FontName','Times New Roma')
xlabel('time/s','FontName','Times New Roma')
title(titlename, 'FontSize', 10, 'FontWeight', 'normal')

% set time0name
xtklabels = xticklabels();
xtklabels{cellfun(@(x) strcmp(x, '0'), xtklabels)} = time0name;
xticklabels(xtklabels);

% time0name line
plot([0 0], ylim, '--r', 'LineWidth',2)

colorbar