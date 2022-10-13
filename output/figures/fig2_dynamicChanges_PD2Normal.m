function fig2_dynamicChanges_PD2Normal(varargin)
%
%   Usage:
%       fig2_dynamicChanges_PD2Normal('newsavefolder', true)
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


title_prefix_original = '-original-ciCohChanges-';
title_prefix_sig = '-sig-ciCohChanges-';


%% Code start here
disp(['running ' funcname]);
time0name = 'reachonset';
t_AOI = [-1 1];
cbarStr = 'ciCohChange';
cbarTicks = [-1 0  1];
histClim = [-1 1];
for ai = 1 : length(animals)
    animal = animals{ai};
    
    if strcmpi(animal, 'Jo')
        input_folder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm4_uNHP_ciCohChanges_PD2Normal_dynamic');
    end
    if strcmpi(animal, 'Kitty')
        input_folder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'longerTrials','m4_uNHP_ciCohChanges_PD2Normal_dynamic');
    end
    
    ciCohChangesfile = fullfile(input_folder, [animal '-ciCohChanges_PD2normal.mat']);
    load(ciCohChangesfile, 'ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea', 't_selected')
    
    % extract chnsNames
    nchns = size(ciCohChanges, 1);
    chnsNames = T_chnsarea.brainarea;
    chnsNames{cellfun(@(x) contains(x, 'stn'), chnsNames)} = 'STN';
    chnsNames{cellfun(@(x) contains(x, 'gp'), chnsNames)} = 'GP';
    
    % extract data from t_AOI
    mask_tAOI = (t_selected >= t_AOI(1) & t_selected <= t_AOI(2));
    ciCohChanges_tAOI = ciCohChanges(:,:, :, mask_tAOI);
    t_selected_AOI = t_selected(mask_tAOI);
    psedociCohChanges_cell = struct2cell(psedociCohChanges);
    psedociCohChanges_AOI = psedociCohChanges_cell(mask_tAOI);
    clear mask_tAOI psedociCohChanges_cell
    
    
   
    % show tags
    show_ylabel = true;
    show_xlabel = true;
    show_titlename = false;
    show_colorbar = true;
    if ai ~= 1
        show_ylabel = false;
    end
    if ai ~= length(animals)
        show_colorbar = false;
    end
    
    
    %%% plot original dynamic ciCohChanges
    title_prefix = [animal title_prefix_original];
    for chi = 1 : nchns -1
        chnnamei = chnsNames{chi};
        for chj = chi + 1 : nchns
            chnnamej = chnsNames{chj};
            ciCohChange_1pair = squeeze(ciCohChanges_tAOI(chi, chj, :, :));
            
            titlename = [title_prefix chnnamei '-' chnnamej];
            
            %%% plot
            ifig = figure('Position', pos_ifig);
            plot_dynamic_ciCohHist_1pair(ciCohChange_1pair, f_selected, t_selected_AOI, time0name, ...
                'fig', ifig, ...
                'histClim', histClim, 'cbarStr', cbarStr, 'cbarTicks', cbarTicks, 'titlename', titlename, ...
                'show_ylabel', show_ylabel, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename, 'show_colorbar', show_colorbar, ...
                'codesavefolder', savecodefolder);
            
            
            % save
            savefilename = titlename;
            print(ifig, fullfile(savefolder, savefilename), '-dpng', '-r1000')
            if ~isempty(copy2folder)
                print(ifig, fullfile(copy2folder, savefilename), '-painters', '-depsc')
            end
            close(ifig)
            
            clear chnnamej ciCohChange_1pair titlename savefilename
        end
    end
    clear

    %%% plot sig dynamic ciCohChanges
    title_prefix = [animal title_prefix_sig];
    
    % extract sigCicohChange_tAOI
    sigCicohChange_tAOI = [];
    for ti = 1 : length(psedociCohChanges_AOI)
        psedociCohChanges_1t = psedociCohChanges_AOI{ti};
        ciCohChange_tAOI_1t = squeeze(ciCohChanges_tAOI(:, :, :, ti));
        sigCicohChange_1t = sigciCoh_extract(psedociCohChanges_1t, ciCohChange_tAOI_1t, 'codesavefolder', savecodefolder);
        
        sigCicohChange_tAOI = cat(4, sigCicohChange_tAOI, sigCicohChange_1t);
        clear tname psedociCohChanges_1t  ciCohChange_tAOI_1t sigCicohChange_1t
    end
    
    for chi = 1 : nchns -1
        chnnamei = chnsNames{chi};
        for chj = chi + 1 : nchns
            chnnamej = chnsNames{chj};
            sigciCohChange_1pair = squeeze(sigCicohChange_tAOI(chi, chj, :, :));
            
            titlename = [title_prefix chnnamei '-' chnnamej '-' pdcond];
            
            % plot
            ifig = figure('Position', pos_ifig);
            plot_dynamic_ciCohHist_1pair(sigciCohChange_1pair, f_selected, t_selected_AOI, time0name, ...
                'fig', ifig, ...
                'histClim', histClim, 'cbarStr', cbarStr, 'cbarTicks', cbarTicks, 'titlename', titlename, ...
                'show_ylabel', show_ylabel, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename, 'show_colorbar', show_colorbar, ...
                'codesavefolder', savecodefolder);
            
            
            % save
            savefilename = titlename;
            print(ifig, fullfile(savefolder, savefilename), '-dpng', '-r1000')
            if ~isempty(copy2folder)
                print(ifig, fullfile(copy2folder, savefilename), '-painters', '-depsc')
            end
            close(ifig)
            
            clear chnnamej ciCohChange_1pair titlename savefilename
        end
    end
    clear title_prefix sigCicohChange_tAOI
    
    
    clear animal input_folder 
    clear ciCohChangesfile
    clear('ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea', 't_selected');
    clear ciCohChanges_tAOI  t_selected_AOI
end






