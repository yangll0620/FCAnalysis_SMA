function fig2a_ccAmpChanges_compCond(varargin)
%
%   Usage:
%       fig2a_ccAmpChanges_compCond('newsavefolder', true)
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
animals = {'Jo'; 'Kitty'};

[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
[~, funcname, ~]= fileparts(codefilepath);

savefolder = fullfile(outputfolder, 'results', 'figures', 'current', funcname);
if(exist(savefolder, 'dir') && newsavefolder)
    disp(['rmdir old savefolder']);
    rmdir(savefolder, 's')
end
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', 'Current',funcname);
if(exist(aisavefolder, 'dir') && newsavefolder)
    disp(['rmdir old aisavefolder']);
    rmdir(aisavefolder, 's')
end
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end
copy2folder = aisavefolder;


savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, "dir")
    disp(['rmdir old savecodefolder']);
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;

%% Code start here
disp(['running ' funcname]);

for ai = 1 : length(animals)
    animal = animals{ai};
    input_folder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm4_uNHP_ccAmpChanges_compCond');
    
    % show_yticklabels and show_colorbar
    show_yticklabels = false;
    show_colorbar = false;
    if(ai == length(animals))
        show_colorbar = true;
    end

    % plot
    plot_ccAmpChanges(pos_ifig, savefilename, animal, input_folder, savefolder, savecodefolder, copy2folder, show_yticklabels, show_colorbar);
    
    clear animal input_folder show_yticklabels show_colorbar
end

function plot_ccAmpChanges(pos_ifig, savefilename, animal, input_folder, savefolder, savecodefolder, copy2folder, show_yticklabels, show_colorbar)
ePhases = {'preMove'; 'earlyReach';  'PeakV'; 'lateReach'};


[tbl_compConds, ~]= compCondEvents_extract(animal);
for compCi = 1: height(tbl_compConds)
    basepd = tbl_compConds.baseCond{compCi};
    comppd = tbl_compConds.compCond{compCi};
    
    for ei = 1 : length(ePhases)
        event = ePhases{ei};        
        ccAmpChangesfile = fullfile(input_folder, [animal '-ccAmpChanges_' comppd '2' basepd '_' event '.mat']);
        
        if ~exist(ccAmpChangesfile, 'file')
            clear event ccAmpChangesfile
            disp([ccAmpChangefile ' not exist']);
            continue;
        end
        
        load(ccAmpChangesfile, 'ccAmpChanges', 'psedoccAmpChanges', 'ccAmp_base', 'ccAmp_comp', 'f_selected', 'T_chnsarea')
        
        % sig
        [sigccAmpChanges]= sigciCoh_extract(psedoccAmpChanges, ccAmpChanges);
        
        
        % remove sig changes where both original ccAmp not sig
        masks_BothNosigs = (ccAmp_base == 0) & (ccAmp_comp == 0);
        sigccAmpChanges(masks_BothNosigs)= 0;
        clear masks_BothNosigs
        
        % flatten
        [sigccAmpChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigccAmpChanges, T_chnsarea);
        [chnPairNames]= chnPairNames_wonum(chnPairNames);
        
        show_titlename = false;
        show_xlabel = false;
        show_xticklabels = false;
        if ei == 1
            show_titlename = true;
        end
        if ei == length(ePhases)
            show_xlabel = true;
            show_xticklabels = true;
        end
        
        % plot subfigure
        plot1PairHist(animal, event, savefilename, chnPairNames, sigccAmpChanges_flatten, pos_ifig, f_selected, comppd, basepd,...
            savecodefolder, savefolder, copy2folder,...
            show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)
        
        clear event align2name ciCohChangesfile
        clear ccAmpChanges psedoccAmpChanges f_selected  T_chnsarea
        clear sigccAmpChanges sigccAmpChanges_flatten chnPairNames
        clear show_titlename show_xlabel show_xticklabels
    end
    
    clear basepd comppd
end


function plot1PairHist(animal, event, savefilename, chnPairNames, sigciCohChanges_flatten, pos_ifig, f_selected, comppd, basepd, ...
    savecodefolder, savefolder, copy2folder,...
    show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)

% plot subfigure
for cpi = 1 : length(chnPairNames)
    chnPairName = chnPairNames{cpi};
    sigciCohChanges_flatten_1pair = sigciCohChanges_flatten(cpi, :);

    ifig = figure('Position', pos_ifig);
    set(ifig, 'PaperUnits', 'points');
    plot_ciCohHistogram_1pair(sigciCohChanges_flatten_1pair, chnPairName, f_selected, [comppd '-' basepd], 'histClim', [-1 1],...
        'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
        'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
        'fig', ifig);

    subsavefolder = fullfile(savefolder, chnPairName);
    if(~exist(subsavefolder,'dir'))
        mkdir(subsavefolder);
    end

    subfilename = [savefilename '-' animal '-' chnPairName comppd '2' basepd '-' event];
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


