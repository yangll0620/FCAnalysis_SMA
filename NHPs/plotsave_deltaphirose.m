function plotsave_deltaphirose(deltaphis_flatten, ciCoh_flatten, chnPairNames, f_selected, titlename_prefix, subtitlename, savefolder, savefile_prefix, savefile_suffix, image_type, varargin)
%   Inputs:
%       
%       Name-Value: 
%           'codesavefolder' - code saved folder         
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

nbins = 10;

fig_left = 50;
fig_bottom = 50;
fig_width = 800;
fig_height = 800;


nchnPairs = length(chnPairNames);
nf = length(f_selected);
for chnPairi = 1 : nchnPairs
    chnPairName = chnPairNames{chnPairi};
    
    subchnpairsavefolder = fullfile(savefolder, chnPairName);
    if ~exist(subchnpairsavefolder, 'dir')
        mkdir(subchnpairsavefolder);
    end
    
    for nfi = 1 : nf
        deltaphi = squeeze(deltaphis_flatten(chnPairi, nfi, :));
        f = f_selected(nfi);
        
        figure;
        set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
        set(gca, 'Position', [0.05 0.05 0.85 0.85])
        polarhistogram(deltaphi, nbins, 'Normalization', 'probability');
        
        titlename = [titlename_prefix ' Trial Phase Diff of ' chnPairName ' at ' num2str(round(f)) 'Hz'];
        title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
        
        % subtitle
        annotation(gcf,'textbox',[0.5 0.017 0.5 0.032], 'String',{subtitlename}, 'LineStyle','none', 'FitBoxToText','off');
        
        % plot icoh if sig
        sig = false;
        icoh = ciCoh_flatten(chnPairi, nfi);
        if(icoh > 0)
            sig = true;
        end
        if sig
            text(0.8, 0.9, ['icoh = ' num2str(round(icoh, 3)) '*'], 'FontSize', 12);
        end
        
        % save
        sigstr = '';
        if sig
            sigstr = ['_sig'];
        end
        savefile =  fullfile(subchnpairsavefolder, [savefile_prefix '_pair' chnPairName '_' num2str(round(f))  'Hz_' savefile_suffix sigstr '.' image_type]);
        saveas(gcf,savefile, image_type);
        
        clear icoh
        close all
    end
end