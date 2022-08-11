function fig2_imCohChanges_compCond(varargin)
%
%   Usage:
%       fig2_imCohChanges_compCond('newsavefolder', true)
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

[~, ~, ~, outputfolder] = exp_subfolders();

[~, funcname, ~]= fileparts(codefilepath);



savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if(exist(savefolder, 'dir') && newsavefolder)
    rmdir(savefolder, 's')
end
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', 'Current',funcname);
if(exist(aisavefolder, 'dir') && newsavefolder)
    rmdir(aisavefolder, 's')
end
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end
copy2folder = aisavefolder;



savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;



%% Code start here
disp(['running ' funcname]);


animals = {'Jo'; 'Kitty'};


for ai = 1 : length(animals)

    animal = animals{ai};

    show_yticklabels = false;
    show_colorbar = false;
    if(ai == length(animals))
        show_colorbar = true;
    end

    if strcmpi(animal, 'Jo')
        plot_JoChanges(pos_ifig, savefilename, savefolder, savecodefolder, copy2folder, show_yticklabels, show_colorbar)

    elseif strcmp(animal, 'Kitty')
        plot_KittyChanges(pos_ifig, savefilename, savefolder, savecodefolder, copy2folder, show_yticklabels, show_colorbar)
    end

    clear animal

end

function plot_KittyChanges(pos_ifig, savefilename, savefolder, savecodefolder, copy2folder, show_yticklabels, show_colorbar)
ePhases = {'preMove'; 'earlyReach';  'PeakV'; 'lateReach'};
animal = 'Kitty';
[~, ~, pipelinefolder, ~] = exp_subfolders();
input_folder = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compCond');

for ei = 1 : length(ePhases)

    % extract sigciCohChanges_flatten
    event = ePhases{ei};
    [~, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, 'pdcond', 'moderate', 'codesavefolder', savecodefolder);

    ciCohChangesfile = fullfile(input_folder, [animal '-ciCohChanges'  '_bnormal' '--moderate' '_' event '_align2' align2name '.mat']);


    if ~exist(ciCohChangesfile, 'file')
        clear event align2name ciCohChangesfile
        continue;
    end

    load(ciCohChangesfile, 'ciCoh_base','ciCoh_comp','psedociCohs_comp','psedociCohs_base','ciCohChanges', 'psedoiCohChanges', 'f_selected', 'T_chnsarea')

    % sig
    [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);


    % remove sig changes where both original ciCoh not sig
    [sigciCoh_base]= sigciCoh_extract(psedociCohs_base, ciCoh_base);
    [sigciCoh_comp]= sigciCoh_extract(psedociCohs_comp, ciCoh_comp);
    masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
    sigciCohChanges(masks_BothNosigs)= 0;
    clear masks_BothNosigs

    % flatten
    [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
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
    plot1PairHist(animal, event, savefilename, chnPairNames, sigciCohChanges_flatten, pos_ifig, f_selected, savecodefolder, savefolder, copy2folder,...
    show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)

    clear event align2name ciCohChangesfile
    clear ciCohChanges psedoiCohChanges f_selected  T_chnsarea
    clear sigciCohChanges sigciCohChanges_flatten chnPairNames
    clear show_titlename show_xlabel show_xticklabels
end



function plot_JoChanges(pos_ifig, savefilename, savefolder, savecodefolder, copy2folder, show_yticklabels, show_colorbar)
ePhases = {'preMove'; 'earlyReach';  'PeakV'; 'lateReach'};
animal = 'Jo';
[~, ~, pipelinefolder, ~] = exp_subfolders();
input_folder = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm5_imCohChanges_PD2normal');

for ei = 1 : length(ePhases)

    % extract sigciCohChanges_flatten
    event = ePhases{ei};
    [~, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(event, 'Jo', 'codesavefolder', savecodefolder);

    ciCohChangesfile = fullfile(input_folder, ['Jo-ciCohChanges_PD2Normal__PD2Normal_' event '_align2' align2name '.mat']);

    if ~exist(ciCohChangesfile, 'file')
        clear event align2name ciCohChangesfile
        continue;
    end

    load(ciCohChangesfile, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'f_selected', 'T_chnsarea')

    % sig
    [sigciCohChanges]= sigciCoh_extract(psedociCohChanges, ciCohChanges);


    % remove sig changes where both original ciCoh not sig
    [sigciCoh_base]= sigciCoh_extract(psedociCohs.normal, ciCohs.normal);
    [sigciCoh_comp]= sigciCoh_extract(psedociCohs.PD, ciCohs.PD);
    masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
    sigciCohChanges(masks_BothNosigs)= 0;
    clear masks_BothNosigs

    % flatten
    [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
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

    plot1PairHist(animal, event, savefilename, chnPairNames, sigciCohChanges_flatten, pos_ifig, f_selected, savecodefolder, savefolder, copy2folder,...
    show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)



    clear event align2name ciCohChangesfile
    clear ciCohChanges psedoiCohChanges f_selected  T_chnsarea
    clear sigciCohChanges sigciCohChanges_flatten chnPairNames
    clear show_titlename show_xlabel show_xticklabels
end


function plot1PairHist(animal, event, savefilename, chnPairNames, sigciCohChanges_flatten, pos_ifig, f_selected, savecodefolder, savefolder, copy2folder,...
    show_xticklabels, show_yticklabels, show_xlabel, show_titlename, show_colorbar)

% plot subfigure
for cpi = 1 : length(chnPairNames)
    chnPairName = chnPairNames{cpi};
    sigciCohChanges_flatten_1pair = sigciCohChanges_flatten(cpi, :);

    ifig = figure('Position', pos_ifig);
    set(ifig, 'PaperUnits', 'points');
    plot_ciCohHistogram_1pair(sigciCohChanges_flatten_1pair, chnPairName, f_selected, 'PD-Normal', 'histClim', [-1 1],...
        'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
        'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
        'fig', ifig);

    subsavefolder = fullfile(savefolder, chnPairName);
    if(~exist(subsavefolder,'dir'))
        mkdir(subsavefolder);
    end

    subfilename = [savefilename '-' animal '-' chnPairName '-PD2Normal' '-' event];
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


