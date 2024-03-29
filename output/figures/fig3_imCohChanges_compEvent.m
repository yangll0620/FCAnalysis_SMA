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
addParameter(p, 'pos_ifig', [50 50 450 150], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
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

input_folder_J_Normal = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');
input_folder_J_PD = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm5_imCohChanges_PD_compEvent');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compEvent');



savefolder = fullfile(outputfolder, 'results', 'figures', 'current', funcname);

% copy to ai savefolder
aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', 'Current',funcname);
if(exist(aisavefolder, 'dir') && newsavefolder)
    rmdir(aisavefolder, 's')
end
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end
copy2folder = aisavefolder;


% save code
savecodefolder = fullfile(savefolder, 'code');
if(exist(savefolder, 'dir') && newsavefolder)
    rmdir(savefolder, 's')
end
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;




%% Code start here
disp(['running ' funcname]);

baseEvent = 'preMove';

ePhases_J = {'earlyReach'; };
ePhases_K = {'earlyReach'; };

cond_J = {'Normal'; 'PD'};
cond_K = {'Normal'; 'Moderate'};


animals = {'Jo'; 'Kitty'};

close all


for ai = 1 : length(animals)
    animal = animals{ai};

    if strcmpi(animal, 'Jo')
        pdconds = cond_J;
        ePhases = ePhases_J;
        
    elseif strcmp(animal, 'Kitty')
        input_folder = input_folder_K;
        pdconds = cond_K;
        ePhases = ePhases_K; 
    end

    for ci = 1 : length(pdconds)
        pdcond = pdconds{ci};
        
        if strcmpi(animal, 'Jo') && strcmpi(pdcond, 'Normal')
            input_folder = input_folder_J_Normal;
        end
        if strcmpi(animal, 'Jo') && strcmpi(pdcond, 'PD')
            input_folder = input_folder_J_PD;
        end
        if strcmpi(animal, 'Kitty') 
            input_folder = input_folder_K;
        end

        show_yticklabels = false;
        show_colorbar = false;
        if ci == length(pdconds) && ai == length(animals)
            show_colorbar = true;
        end


        for ei = 1 : length(ePhases)
            
            % extract sigciCohChanges_flatten
            compEvent = ePhases{ei};
            
            if strcmpi(animal, 'Jo') && strcmpi(pdcond, 'PD')
                ciCohChangesfile = fullfile(input_folder, ['Jo-ciCohChanges_PD__PD_b' baseEvent '--' compEvent '.mat']);
                load(ciCohChangesfile, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'f_selected', 'T_chnsarea');
                
                ciCoh_base = ciCohs.(baseEvent);
                ciCoh_comp = ciCohs.(compEvent);
                psedociCohs_base = psedociCohs.(baseEvent);
                psedociCohs_comp = psedociCohs.(compEvent);
                psedoiCohChanges = psedociCohChanges;
                
                clear ciCohs psedociCohs psedociCohChanges
            else
                [~, ~, align2name_comp] = SKT_EventPhase_align2_tAOI_extract(compEvent, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
                ciCohChangesfile = fullfile(input_folder, [animal '-ciCohChanges_' pdcond '_b' baseEvent '--' compEvent '_align2' align2name_comp '.mat']);
                
                
                load(ciCohChangesfile, 'ciCoh_base','ciCoh_comp','psedociCohs_comp','psedociCohs_base','ciCohChanges', 'psedoiCohChanges', 'f_selected', 'T_chnsarea');
                clear align2name_comp
                
            end
            
            

            [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);


            % remove sig changes where both original ciCoh not sig
            [sigciCoh_base]= sigciCoh_extract(psedociCohs_base, ciCoh_base);
            [sigciCoh_comp]= sigciCoh_extract(psedociCohs_comp, ciCoh_comp);
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
        
        clear input_folder show_yticklabels show_colorbar
    end

    clear pdconds ePhases
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
