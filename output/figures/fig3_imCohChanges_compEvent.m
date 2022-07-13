function fig3_imCohChanges_compEvent(varargin)
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

input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');


savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
    
end

% copy to ai savefolder
aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', funcname);
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

baseEvent = 'preMove';

ePhases_J = {'earlyReach';  'PeakV'; 'lateReach'};
ePhases_K = {'earlyReach';  'PeakV'; 'lateReach'};

cond_J = {'Normal'; 'Mild';  'Moderate'};
cond_K = {'Normal'; 'Moderate'};


animals = {'Jo'; 'Kitty'};

close all

histClim_changes = [-1 1];
cbarStr_changes = 'ciCohChange';
cbarTicks_changes = [-1 0 1];

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
        if ci == 1 && ai == 1
            show_yticklabels = true;
        end
        if ci == length(pdconds) && ai == length(animals)
            show_colorbar = true;
        end


        for ei = 1 : length(ePhases)
            
            % extract sigciCohChanges_flatten
            compEvent = ePhases{ei};
    
            [~, ~, align2name_comp] = SKT_EventPhase_align2_tAOI_extract(compEvent, animal, pdcond, 'codesavefolder', savecodefolder);
            
            ciCohChangesfile = fullfile(input_folder, [animal '-ciCohChanges_' pdcond '_b' baseEvent '--' compEvent '_align2' align2name_comp '.mat']);
            load(ciCohChangesfile, 'ciCoh_base','ciCoh_comp','psedociCohs_comp','psedociCohs_base','ciCohChanges', 'psedoiCohChanges', 'f_selected', 'T_chnsarea');

            
            [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);


            % remove sig changes where both original ciCoh not sig
            [sigciCoh_base]= sigciCoh_extract(psedociCohs_base, ciCoh_base);
            [sigciCoh_comp]= sigciCoh_extract(psedociCohs_comp, ciCoh_comp);
            masks_BothNosigs = (sigciCoh_base == 0) & (sigciCoh_comp == 0);
            sigciCohChanges(masks_BothNosigs)= 0;
            clear masks_BothNosigs sigciCoh_base sigciCoh_comp

            [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);

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
            ifig = figure('Position', pos_ifig);
            set(ifig, 'PaperUnits', 'points');
            plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, pdcond, 'histClim', histClim_changes,...
                'codesavefolder', savecodefolder, 'cbarStr', cbarStr_changes, 'cbarTicks', cbarTicks_changes, ...
                'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
                'fig', ifig);

            
            subfilename = [savefilename '-' animal '-' compEvent '-B' baseEvent '-' pdcond]; 
            print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
            print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')

            if ~isempty(copy2folder)
                print(ifig, fullfile(copy2folder, subfilename), '-painters', '-depsc')
            end

            close(ifig)


            clear compEvent align2name_comp
            clear('ciCoh_base','ciCoh_comp','psedociCohs_comp','psedociCohs_base','ciCohChanges', 'psedoiCohChanges', 'f_selected', 'T_chnsarea'); 
            clear sigciCohChanges sigciCohChanges_flatten chnPairNames
        end
    end

    clear input_folder pdconds ePhases
end

