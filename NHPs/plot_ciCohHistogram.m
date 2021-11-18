function plot_ciCohHistogram(ciCoh_flatten, chnPairNames, f_selected, titlename, varargin)
%
%   Inputs:
%       
%
%

if nargin < 1
    histClim = [0 1];
else
    histClim = varargin{1};
end

fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;


% plot
figure;
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
imagesc(ciCoh_flatten)
colormap(jet)
set(gca, 'Position', [0.09 0.05 0.9 0.88])
[npairs, nf] = size(ciCoh_flatten);
xticks([1:nf])
xticklabels(round(f_selected,2))
yticks([1:npairs]);
set(gca,'YTickLabel',chnPairNames,'fontsize',12,'FontWeight','bold')
xlabel('freqs')
title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
set(gca,'CLim', histClim)
colorbar



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