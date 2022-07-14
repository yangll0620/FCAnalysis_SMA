function fig4_imCohChanges_FastBasedSlowReach(varargin)
%   
%   Usage:
%       fig4_imCohChanges_FastBasedSlowReach('pos_ifig', [50 50 400 200])
%
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the reachtime statistical figure [left bottom fig_width fig_height], default [50 50 400 200]
%  


% parse params
p = inputParser;
addParameter(p, 'pos_ifig', [50 50 400 200], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));


parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;



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

animals = {'Kitty'};

ePhases = {'earlyReach';  'PeakV'; 'lateReach'};
for ai = 1 : length(animals)

    animal = animals{ai};
    input_folder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm4_uNHP_imCohChanges_FastbasedSlowReach');

    ciCoh_Changes_file = fullfile(input_folder, 'ciCohs_FastbasedSlowReach_ciCohChanges.mat');

    load(ciCoh_Changes_file, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea');

    show_yticklabels = false;
    show_colorbar = true;
    if ai ==1 
        show_yticklabels = true;
    end
    if ai == length(animals)
        show_colorbar = true;
    end

    for ei = 1 : length(ePhases)
        ePhase = ePhases{ei};

        ciCoh_slow = ciCohs.slowReach.(ePhase);
        ciCoh_fast = ciCohs.fastReach.(ePhase);
        ciCohChange_fastBslow= ciCohChanges.(ePhase);

        psedoCicoh_slow = psedociCohs.(ePhase).slowReach;
        psedoCicoh_fast = psedociCohs.(ePhase).fastReach;
        psedociCohChange_fastBslow = psedociCohChanges.(ePhase);


        % sigciCohChanges
        [sigciCohChanges]= sigciCoh_extract(psedociCohChange_fastBslow, ciCohChange_fastBslow);


        % remove sig changes where both original ciCoh not sig
        ciCoh_base = ciCoh_slow;
        ciCoh_comp = ciCoh_fast;
        psedociCohs_base = psedoCicoh_slow;
        psedociCohs_comp = psedoCicoh_fast;
        [sigciCoh_base]= sigciCoh_extract(psedociCohs_base, ciCoh_base);
        [sigciCoh_comp]= sigciCoh_extract(psedociCohs_comp, ciCoh_comp);
        masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
        sigciCohChanges(masks_BothNosigs)= 0;
        clear masks_BothNosigs 
        clear ciCoh_base ciCoh_comp psedociCohs_base psedociCohs_comp



        % flatten
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);


        show_xlabel = false;
        show_xticklabels = false;
        show_titlename = false;
        if ei == length(ePhases)
            show_xlabel = true;
            show_xticklabels = true;
            show_titlename = true;
        end


        % plot subfigure
        ifig = figure('Position', pos_ifig);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, 'Fast-Slow Trials', 'histClim', [-1 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig);

        subfilename = [savefilename '-' animal '-' ePhase]; 
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')

        if ~isempty(copy2folder)
            print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
        end

        close(ifig)


        clear ciCoh_slow ciCoh_fast ciCohChange_fastBslow
        clear psedoCicoh_slow psedoCicoh_fast psedociCohChange_fastBslow
    end
    

end