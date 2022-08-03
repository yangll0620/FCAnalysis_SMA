function fig6_imCohChanges_Reachfreeze2ReachPhases(varargin)
%   
%   Usage:
%       fig6_imCohChanges_Reachfreeze2ReachPhases('pos_ifig', [50 50 400 200], 'plot_ciCohs_ReachFreeze', false, 'plot_ciCohChanges_ReachFreeze2Reach', false, 'plot_ciCohChanges_ReachFreezeAlongTime', false,'plot_ciCoh_Reach', false)
%
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the reachtime statistical figure [left bottom fig_width fig_height], default [50 50 400 200]
%           'plot_ciCohs_ReachFreeze' - tag plotting ciCoh of early/middle/late/after reach freezes (default true)
%           'plot_ciCohChanges_ReachFreeze2Reach' - tag plotting ciCohChanges of early/middle/late/after reach freeze relative to preMove/earlyReach(default false)
%           'plot_ciCohChanges_ReachFreezeAlongTime' - tag plotting ciCohChanges of reach freeze along time, e.g. middleFreeze2earlyFreeze (default true)
%           'plot_ciCoh_Reach' - tag plotting ciCoh of preMove/earlyReach


% parse params
p = inputParser;
addParameter(p, 'pos_ifig', [50 50 400 200], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'plot_ciCohs_ReachFreeze', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plot_ciCohChanges_ReachFreeze2Reach', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plot_ciCohChanges_ReachFreezeAlongTime', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plot_ciCoh_Reach', true, @(x) assert(islogical(x) && isscalar(x)));


parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;
plot_ciCohs_ReachFreeze = p.Results.plot_ciCohs_ReachFreeze;
plot_ciCohChanges_ReachFreeze2Reach = p.Results.plot_ciCohChanges_ReachFreeze2Reach;
plot_ciCohChanges_ReachFreezeAlongTime = p.Results.plot_ciCohChanges_ReachFreezeAlongTime;
plot_ciCoh_Reach = p.Results.plot_ciCoh_Reach;

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

input_folder = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm5_imCohChanges_reachFreeze2ReachPhases');
ciCoh_Changes_file = fullfile(input_folder, 'ciCohs-Changes-reachFreeze2reachPhases.mat');


input_folder2 = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm6_imCohChanges_reachFreeze_alongTime');
ciCoh_Changes_file2 = fullfile(input_folder2, 'ciCohsChanges-reachFreeze-alongTime.mat');

savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', funcname);
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


% plot ciCohs.ReachFreeze.earlyFreeze....
if plot_ciCohs_ReachFreeze
    load(ciCoh_Changes_file, 'ciCohs','psedociCohs', 'f_selected',  'T_chnsarea')
    reachfreezeTypes = fieldnames(ciCohs.ReachFreeze);
    
    for fri = 1 : length(reachfreezeTypes)
        subfreezeType = reachfreezeTypes{fri};

        ciCoh = ciCohs.ReachFreeze.(subfreezeType);
        psedociCoh = psedociCohs.ReachFreeze.(subfreezeType);


        [sigciCoh]= sigciCoh_extract(psedociCoh, ciCoh);
        [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);

        show_xticklabels = false;
        show_xlabel = false;
        show_titlename = false;
        show_yticklabels = false;
        show_colorbar = false;
        if fri == 1
            show_yticklabels = true;
        end
        if fri == length(reachfreezeTypes)
            show_colorbar = true;
        end


        % plot and save ciCoh Histogram image
        ifig = figure('Position', pos_ifig);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, ['reachFreeze ' subfreezeType], 'histClim', [0 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCoh', 'cbarTicks', [0 0.5 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig);

        subfilename = [savefilename '-' subfreezeType];
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')

        if ~isempty(copy2folder)
            print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
        end

        close(ifig)


        clear subfreezeType
        clear ciCoh psedociCoh sigciCoh sigciCoh_flatten
        clear show_xticklabels  show_xlabel
        clear show_titlename show_yticklabels show_colorbar
        clear ifig
    end

    clear('ciCohs','psedociCohs', 'f_selected',  'T_chnsarea');
    clear reachfreezeTypes
end


% plot ciCohChanges_bpreMove/bearlyReach
if plot_ciCohChanges_ReachFreeze2Reach
    load(ciCoh_Changes_file, 'ciCohs','psedociCohs','ciCohChanges','psedociCohChanges', 'f_selected',  'T_chnsarea')
    eBasePhases = fieldnames(ciCohChanges);
    reachfreezeTypes = fieldnames(ciCohs.ReachFreeze);

    for ebi = 1 : length(eBasePhases)
        eBasePhase = eBasePhases{ebi};

        ciCoh_base = ciCohs.(eBasePhase(3:end));
        psedociCohs_base = psedociCohs.Reach.(eBasePhase(3:end));
        [sigciCoh_base]= sigciCoh_extract(psedociCohs_base, ciCoh_base);

        show_xticklabels = false;
        show_xlabel = false;
        if ebi == length(eBasePhases)
            show_xticklabels = true;
            show_xlabel = true;
        end

        for fri = 1 : length(reachfreezeTypes)
            subfreezeType = reachfreezeTypes{fri};

            ciCoh_comp = ciCohs.ReachFreeze.(subfreezeType);
            psedociCohs_comp = psedociCohs.ReachFreeze.(subfreezeType);


            ciCohChange = ciCohChanges.(eBasePhase).ReachFreeze.(subfreezeType);
            psedociCohChange = psedociCohChanges.(['b' eBasePhase(3:end)]).ReachFreeze.(subfreezeType);
            [sigciCohChanges]= sigciCoh_extract(psedociCohChange, ciCohChange);


            % remove sig changes where both original ciCoh not sig
            [sigciCoh_comp]= sigciCoh_extract(psedociCohs_comp, ciCoh_comp);
            masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
            sigciCohChanges(masks_BothNosigs)= 0;



            [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);

            show_titlename = false;
            show_yticklabels = false;
            show_colorbar = false;
            if fri == 1
                show_yticklabels = true;
            end
            if fri == length(reachfreezeTypes)
                show_colorbar = true;
            end


            % plot and save ciCoh Histogram image
            ifig = figure('Position', pos_ifig);
            set(ifig, 'PaperUnits', 'points');
            plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, ['reachFreeze ' subfreezeType '-' eBasePhase], 'histClim', [-1 1],...
                'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
                'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
                'fig', ifig);

            subfilename = [savefilename '-' eBasePhase '-' subfreezeType];
            print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
            print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')

            if ~isempty(copy2folder)
                print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
            end

            close(ifig)

            clear subfreezeType ciCoh_comp psedociCohs_comp
            clear ciCohChange psedociCohChange sigciCohChanges
            clear sigciCoh_comp masks_BothNosigs
            clear show_titlename show_yticklabels show_colorbar
            clear ifig subfilename
        end

        clear eBasePhase ciCoh_base psedociCohs_base sigciCoh_base
        clear show_xticklabels show_xlabel
    end

    clear('ciCohs','psedociCohs','ciCohChanges','psedociCohChanges', 'f_selected',  'T_chnsarea');
    clear eBasePhases reachfreezeTypes
end



% plot ciCohChanges_middle2EarlyFreeze
if plot_ciCohChanges_ReachFreezeAlongTime
    load(ciCoh_Changes_file2, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'f_selected', 'T_chnsarea');
    reachfreezeTypes = fieldnames(ciCohs.ReachFreeze);

    for fri_base = 1 : length(reachfreezeTypes)-1
        subfreezeType_base = reachfreezeTypes{fri_base};

        ciCoh_base = ciCohs.ReachFreeze.(subfreezeType_base);
        psedociCohs_base = psedociCohs.ReachFreeze.(subfreezeType_base);

        [sigciCoh_base]= sigciCoh_extract(psedociCohs_base, ciCoh_base);

        show_xticklabels = true;
        show_xlabel = true;

        for fri_comp = fri_base + 1 : length(reachfreezeTypes)
            subfreezeType_comp = reachfreezeTypes{fri_comp};

            ciCoh_comp = ciCohs.ReachFreeze.(subfreezeType_comp);
            psedociCohs_comp = psedociCohs.ReachFreeze.(subfreezeType_comp);

            ciCohChange = ciCohChanges.ReachFreeze.([subfreezeType_comp '2' subfreezeType_base]);
            psedociCohChange = psedociCohChanges.ReachFreeze.([subfreezeType_comp '2' subfreezeType_base]);


            % sig
            [sigciCohChanges]= sigciCoh_extract(psedociCohChange, ciCohChange);

            % remove sig changes where both original ciCoh not sig
            [sigciCoh_comp]= sigciCoh_extract(psedociCohs_comp, ciCoh_comp);
            masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
            sigciCohChanges(masks_BothNosigs)= 0;

            %flatten
            [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);


            show_titlename = true;
            show_yticklabels = false;
            show_colorbar = false;
            if fri_comp == 2
                show_yticklabels = true;
            end
            if fri_comp == length(reachfreezeTypes)
                show_colorbar = true;
            end


            % plot and save ciCoh Histogram image
            ifig = figure('Position', pos_ifig);
            set(ifig, 'PaperUnits', 'points');
            plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, [subfreezeType_comp '-' subfreezeType_base], 'histClim', [-1 1],...
                'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
                'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
                'fig', ifig);

            subfilename = [savefilename '-' subfreezeType_comp '2' subfreezeType_base];
            print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
            print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')

            if ~isempty(copy2folder)
                print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
            end

            close(ifig)

            clear subfreezeType_comp ciCoh_comp psedociCohs_comp
            clear ciCohChange psedociCohChange sigciCohChanges sigciCohChanges_flatten chnPairNames
            clear sigciCoh_comp masks_BothNosigs
            clear show_titlename show_yticklabels show_colorbar
            clear ifig subfilename
        end

        clear subfreezeType_base ciCoh_base psedociCohs_base sigciCoh_base
        clear show_xticklabels show_xlabel
    end
    
    clear fri_base
    clear('ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'T_chnsarea');
    clear reachfreezeTypes
end



% plot ciCoh_preMove/earlyReach
if plot_ciCoh_Reach
    load(ciCoh_Changes_file2, 'ciCohs', 'psedociCohs', 'f_selected', 'T_chnsarea');
    reachPhases = {'preMove', 'earlyReach'};

    for ri = 1 : length(reachPhases)
        reachPhase = reachPhases{ri};
    
        ciCoh = ciCohs.(reachPhase);
        psedociCoh = psedociCohs.Reach.(reachPhase);  
    
    
        [sigciCoh]= sigciCoh_extract(psedociCoh, ciCoh);
        [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);
    
    
    
        show_xticklabels = true;
        show_xlabel = true;
        show_titlename = true;
        show_yticklabels = false;
        show_colorbar = false;
        if ri == 1
            show_yticklabels = true;
        end
        if ri == length(reachPhases)
            show_colorbar = true;
        end
    
    
    
        % plot and save 
        ifig = figure('Position', pos_ifig);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, [reachPhase], 'histClim', [0 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCoh', 'cbarTicks', [0 0.5 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig);
    
        subfilename = [savefilename '-' reachPhase];
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
    
        if ~isempty(copy2folder)
            print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
        end
    
        close(ifig)
    
    
    
        clear reachPhase ciCoh psedociCoh 
        clear sigciCoh  sigciCoh_flatten chnPairNames
        clear show_xticklabels  show_xlabel
        clear show_titlename show_yticklabels show_colorbar
        clear ifig subfilename
    end
    clear reachPhases ri
    clear('ciCohs', 'psedociCohs', 'T_chnsarea');
end