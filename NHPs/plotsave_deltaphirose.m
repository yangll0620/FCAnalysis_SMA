function plotsave_deltaphirose(deltaphis_flatten, ciCoh_flatten, chnPairNames, f_selected, titlename_prefix, subtitlename, savefolder, savefile_prefix, savefile_suffix, image_type, varargin)
%   Inputs:
%       
%       Name-Value: 
%           'codesavefolder' - code saved folder
%           'plotNoSig' - plot no sig tag(default false)
%           'roseRLim' - RLim
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'plotNoSig', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'roseRLim', 'auto', @(x) assert(isnumeric(x) && isvector(x) && length(x)==2, 'roseRLim should be a vector with 2 numbers'))
parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
plotNoSig = p.Results.plotNoSig;
if ~isequal(p.Results.roseRLim, 'auto')
    setRoseRLim = true;
    roseRLim = reshape(p.Results.roseRLim, 1, 2);
else
    setRoseRLim = false;
end


% copy code to savefolder if not empty
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
    
    
    [~, maxnfi]= max(ciCoh_flatten(chnPairi,:));
    for nfi = 1 : nf
        icoh = ciCoh_flatten(chnPairi, nfi);
        sig = false;
        if(icoh > 0)
            sig = true;
        end
        if ~plotNoSig && ~sig
            clear icoh sig
            continue;
        end
        
        deltaphi = squeeze(deltaphis_flatten(chnPairi, nfi, :));
        f = f_selected(nfi);
        
        figure;
        set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
        set(gca, 'Position', [0.05 0.05 0.85 0.85])
        polarhistogram(deltaphi, nbins, 'Normalization', 'probability');
        if setRoseRLim
            set(gca, 'RLim', roseRLim)
        end
        
        titlename = [titlename_prefix ' Trial Phase Diff of ' chnPairName ' at ' num2str(round(f)) 'Hz'];
        title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
        
        % subtitle
        annotation(gcf,'textbox',[0.5 0.017 0.5 0.032], 'String',{subtitlename}, 'LineStyle','none', 'FitBoxToText','off');
        
        if sig
            annotation(gcf,'textbox',[0.7 0.8 0.5 0.03], 'String',{['cicoh = ' num2str(round(icoh, 3)) '*']}, 'LineStyle','none', 'FitBoxToText','off');
        end
        
        % save
        sigstr = '';
        if sig
            sigstr = ['_sig'];
        end
        if nfi == maxnfi
            sigstr = [sigstr '_Max'];
        end
        savefile =  fullfile(subchnpairsavefolder, [savefile_prefix '_pair' chnPairName '_' num2str(round(f))  'Hz_' savefile_suffix sigstr '.' image_type]);
        saveas(gcf,savefile, image_type);
        
        clear icoh
        close all
    end
    clear chnPairName subchnpairsavefolder maxnfi
end