function fig3_imCohChanges_JoPD_compEvent(varargin)
%   
%   Usage:
%       fig3_imCohChanges_compEvent('pos_ifig', [50 50 400 200])
%
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the reachtime statistical figure [left bottom fig_width fig_height], default [50 50 400 200]
%  

% parse params
p = inputParser;
addParameter(p, 'pos_ifig', [50 50 410 150], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'newsavefolder', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;
newsavefolder = p.Results.newsavefolder;


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
[~, funcname, ~]= fileparts(codefilepath);

input_folder = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm5_imCohChanges_PD_compEvent');



savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if(exist(savefolder, 'dir') && newsavefolder)
    rmdir(savefolder, 's')
end
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
    
end

% copy to ai savefolder
aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', 'Current', funcname);
if(exist(aisavefolder, 'dir') && newsavefolder)
    rmdir(aisavefolder, 's')
end
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end
copy2folder = aisavefolder;


% save code
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;




%% Code start here
disp(['running ' funcname]);

ePhases = {'earlyReach';  'peakV'; 'lateReach'};
baseEvent = 'preMove';

animals = {'Jo'};

close all

pdcond = 'PD';

for ai = 1 : length(animals)
    animal = animals{ai};


    show_yticklabels = false;
    show_colorbar = false;


    for ei = 1 : length(ePhases)

        % extract sigciCohChanges_flatten
        compEvent = ePhases{ei};


        ciCohChangesfile = fullfile(input_folder, ['Jo-ciCohChanges_PD__PD_b' baseEvent '--' compEvent '.mat']);
        load(ciCohChangesfile, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'f_selected', 'T_chnsarea');


        [sigciCohChanges]= sigciCoh_extract(psedociCohChanges, ciCohChanges);


        % remove sig changes where both original ciCoh not sig
        [sigciCoh_base]= sigciCoh_extract(psedociCohs.(baseEvent), ciCohs.(baseEvent));
        [sigciCoh_comp]= sigciCoh_extract(psedociCohs.(compEvent), ciCohs.(compEvent));
        masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
        sigciCohChanges(masks_BothNosigs)= 0;
        clear masks_BothNosigs sigciCoh_base sigciCoh_comp

        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        [chnPairNames]= chnPairNames_wonum(chnPairNames);

        show_titlename = true;
        show_xlabel = false;
        show_xticklabels = false;
        if ei == length(ePhases)
            show_xlabel = true;
            show_xticklabels = true;
        end

        % plot subfigure
        plotHist_pairbypair(animal, savefilename, chnPairNames, sigciCohChanges_flatten, pos_ifig, f_selected, savecodefolder, savefolder, copy2folder,...
            compEvent, baseEvent, pdcond, ...
            show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)

        clear compEvent align2name_comp
        clear('ciCoh_base','ciCoh_comp','psedociCohs_comp','psedociCohs_base','ciCohChanges', 'psedoiCohChanges', 'f_selected', 'T_chnsarea');
        clear sigciCohChanges sigciCohChanges_flatten chnPairNames
    end

    clear input_folder pdconds ePhases
end

function plotHist_pairbypair(animal, savefilename, chnPairNames, sigciCohChanges_flatten, pos_ifig, f_selected, savecodefolder, savefolder, copy2folder,...
    compEvent, baseEvent, pdcond,...
    show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)

% plot subfigure
for cpi = 1 : length(chnPairNames)
    chnPairName = chnPairNames{cpi};
    sigciCohChanges_flatten_1pair = sigciCohChanges_flatten(cpi, :);

    ifig = figure('Position', pos_ifig);
    set(ifig, 'PaperUnits', 'points');
    plot_ciCohHistogram_1pair(sigciCohChanges_flatten_1pair, chnPairName, f_selected, [compEvent '-' baseEvent], 'histClim', [-1 1],...
        'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
        'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
        'fig', ifig);

    subsavefolder = fullfile(savefolder, chnPairName);
    if(~exist(subsavefolder,'dir'))
        mkdir(subsavefolder);
    end

    subfilename = [savefilename '-' animal '-' chnPairName '-' pdcond '-' compEvent '2' baseEvent];
    print(ifig, fullfile(subsavefolder, subfilename), '-painters', '-depsc')
    print(ifig, fullfile(subsavefolder, subfilename), '-dpng', '-r1000')

    if ~isempty(copy2folder)
        subcopyfolder = fullfile(copy2folder, chnPairName);
        if(~exist(subcopyfolder,'dir'))
            mkdir(subcopyfolder);
        end
        print(ifig, fullfile(subcopyfolder, subfilename), '-painters', '-depsc')
    end

    close(ifig)

    clear chnPairName sigciCohChanges_flatten_1pair
    clear subfilename ifig subcopyfolder subsavefolder
end

