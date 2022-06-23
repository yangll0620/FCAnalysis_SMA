function fig4_imCohChanges_compCond()

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


input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compCond');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm4_fs500Hz_uNHP_imCohChanges_compCond');


savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;



%% Code start here
basepd = 'normal';

compconds_J = cond_cell_extract('Jo');
compconds_K = cond_cell_extract('Kitty');
compconds_J(strcmp(compconds_J, basepd)) = [];
compconds_K(strcmp(compconds_K, basepd)) = [];



ePhases_J = {'preMove'; 'earlyReach';  'lateReach'};
ePhases_K = {'preMove'; 'earlyReach';  'PeakV'; 'lateReach'};

animals = {'Jo'; 'Kitty'};


for ai = 1 : length(animals)

    animal = animals{ai};

    if strcmpi(animal, 'Jo')
        input_folder = input_folder_J;
        ePhases = ePhases_J;
        compconds = compconds_J;
        
    elseif strcmp(animal, 'Kitty')
        input_folder = input_folder_K;
        ePhases = ePhases_K;
        compconds = compconds_K;
    end

    ciCohChangesfile_prefix =[animal '-ciCohChanges'];

    show_yticklabels = false;
    show_colorbar = false;
    nconds = length(compconds);
    for ci = 1 : nconds
        comppd = compconds{ci};

        if ci == 1
            show_yticklabels = true;
        end
        if ci == nconds
            show_colorbar = true;
        end


        show_titlename = false;
        show_xlabel = false;
        show_xticklabels = false;
        nphases = length(ePhases);
        for npi = 1 : nphases
        
            % extract sigciCohChanges_flatten
            event = ePhases{npi};
            [~, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, comppd, 'codesavefolder', savecodefolder);
            ciCohChangesfile = fullfile(input_folder, [ciCohChangesfile_prefix  '_b' basepd '--' comppd '_' event '_align2' align2name '.mat']);
            if ~exist(ciCohChangesfile, 'file')
                clear event align2name ciCohChangesfile
                continue;
            end
            load(ciCohChangesfile, 'ciCohChanges', 'psedoiCohChanges', 'f_selected',  'T_chnsarea')
            [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
            [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        
        
            if npi == 1
               show_titlename = true;
            end
            if npi == nphases
               show_xlabel = true;
               show_xticklabels = true;
            end
        
            % plot subfigure
            ifig = figure('Position', [50 50 400 200]);
            set(ifig, 'PaperUnits', 'points');
            plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, [comppd '-' basepd], 'histClim', [-1 1],...
                'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
                'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
                'fig', ifig);
            
            subfilename = [savefilename '-' animal '-' comppd '-B' basepd '-' event]; % 'Jo-mild-Bnormal-preMove'
            print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
            print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
            close(ifig)
        
          
            clear event align2name ciCohChangesfile
            clear ciCohChanges psedoiCohChanges f_selected  T_chnsarea
            clear sigciCohChanges sigciCohChanges_flatten chnPairNames
        end

        clear comppd
        clear show_titlename show_xlabel show_xticklabels
    end
    clear animal input_folder ePhases compconds nconds
    clear show_yticklabels show_colorbar
    clear ciCohChangesfile_prefix
    
end
