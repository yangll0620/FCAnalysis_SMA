function fig3a_imCohChanges_compEvent(varargin)
%   
%   Usage:
%       fig3a_imCohChanges_compEvent('pos_ifig', [50 50 400 200])
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

input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_uNHP_ccAmpChanges_compEvent');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_uNHP_ccAmpChanges_compEvent');


savefolder = fullfile(outputfolder, 'results', 'figures', 'current', funcname);
if(exist(savefolder, 'dir') && newsavefolder)
    disp(['rmdir old savefolder']);
    rmdir(savefolder, 's')
end
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

% copy to ai savefolder
aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', 'Current',funcname);
if(exist(aisavefolder, 'dir') && newsavefolder)
    disp(['rmdir old aisavefolder']);
    rmdir(aisavefolder, 's')
end
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end
copy2folder = aisavefolder;


% save code
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, "dir")
    disp(['rmdir old savecodefolder']);
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;

%% Code start here
disp(['running ' funcname]);

baseEvent = 'preMove';

ePhases_J = {'earlyReach';  'PeakV'; 'lateReach'};
ePhases_K = {'earlyReach';  'PeakV'; 'lateReach'};

cond_J = {'normal'; 'mild';'moderate'};
cond_K = {'normal'; 'moderate'};


animals = {'Jo'; 'Kitty'};

close all


for ai = 1 : length(animals)
    animal = animals{ai};

    if strcmpi(animal, 'Jo')
        input_folder = input_folder_J;
        pdconds = cond_J;
        ePhases = ePhases_J;
        
    elseif strcmp(animal, 'Kitty')
        input_folder = input_folder_K;
        pdconds = cond_K;
        ePhases = ePhases_K; 
    end

    for ci = 1 : length(pdconds)
        pdcond = pdconds{ci};
        
        show_yticklabels = false;
        show_colorbar = false;
        if ci == length(pdconds) && ai == length(animals)
            show_colorbar = true;
        end
        
        for ei = 1 : length(ePhases)
            compEvent = ePhases{ei};
            ccAmpChangesfile = fullfile(input_folder, [animal '-ccAmpChanges_' pdcond '_' compEvent '2' baseEvent '.mat']);
            
            if ~exist(ccAmpChangesfile, 'file')
                disp([ccAmpChangesfile ' not exist']);
                clear compEvent ccAmpChangesfile
                continue;
            end
            
            load(ccAmpChangesfile, 'ccAmpChanges', 'psedoccAmpChanges', 'ccAmp_base', 'ccAmp_comp', 'f_selected', 'T_chnsarea')

            
            [sigccAmpCohChanges]= sigciCoh_extract(psedoccAmpChanges, ccAmpChanges);


            % remove sig changes where both original ccAmp not sig
            masks_BothNosigs = (ccAmp_base == 0) & (ccAmp_comp == 0);
            sigccAmpCohChanges(masks_BothNosigs)= 0;
            clear masks_BothNosigs 

            % flatten
            [sigccAmpChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigccAmpCohChanges, T_chnsarea);
            [chnPairNames]= chnPairNames_wonum(chnPairNames);

            show_titlename = true;
            show_xlabel = false;
            show_xticklabels = false;
            if ei == length(ePhases)
               show_xlabel = true;
               show_xticklabels = true;
            end
        
            % plot subfigure
            plotHist_pairbypair(animal, savefilename, chnPairNames, sigccAmpChanges_flatten, pos_ifig, f_selected, savecodefolder, savefolder, copy2folder,...
                compEvent, baseEvent, pdcond, ...
                show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)


            clear compEvent 
            clear('ccAmpChanges', 'psedoccAmpChanges', 'ccAmp_base', 'ccAmp_comp', 'f_selected', 'T_chnsarea')
            clear sigccAmpCohChanges sigccAmpChanges_flatten chnPairNames
        end
        clear pdcond show_yticklabels  show_colorbar
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
