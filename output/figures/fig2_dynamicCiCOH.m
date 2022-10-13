function fig2_dynamicCiCOH(varargin)
%
%   Usage:
%       fig2_dynamicCiCOH('newsavefolder', true)
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
cond_cell = {'normal', 'PD'};

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


title_prefix_original = '-original-ciCoh-';
title_prefix_sig = '-sig-ciCoh-';

%% Code start here
disp(['running ' funcname]);
time0name = 'reachonset';
t_AOI = [-1 1];
cbarStr = 'ciCoh';
cbarTicks = [0  0.5 1];
histClim = [0 1];
for ai = 1 : length(animals)
    animal = animals{ai};
    
    if strcmpi(animal, 'Jo')
        input_folder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm3_uNHP_ciCOH_dynamic');
    end
    if strcmpi(animal, 'Kitty')
        input_folder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'longerTrials','m3_uNHP_ciCOH_dynamic');
    end
    
    
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        ciCohfile = fullfile(input_folder, [animal '-ciCoh_Dynamic__' pdcond '.mat']);
        
        % ciCoh: nchns * nchns * nf * nt
        load(ciCohfile, 'ciCoh', 'psedociCohs', 'f_selected',  'T_chnsarea', 't_selected');
        
        % extract chnsNames
        nchns = size(ciCoh, 1);
        chnsNames = T_chnsarea.brainarea;
        chnsNames{cellfun(@(x) contains(x, 'stn'), chnsNames)} = 'STN';
        chnsNames{cellfun(@(x) contains(x, 'gp'), chnsNames)} = 'GP';
        
        % extract data from t_AOI
        mask_tAOI = (t_selected >= t_AOI(1) & t_selected <= t_AOI(2));
        ciCoh_tAOI = ciCoh(:,:, :, mask_tAOI);
        t_selected_AOI = t_selected(mask_tAOI);
        psedociCohs_cell = struct2cell(psedociCohs);
        psedociCohs_AOI = psedociCohs_cell(mask_tAOI);
        clear mask_tAOI psedociCohs_cell
        
        
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
        
        
        %%% plot original dynamic ciCoh
        title_prefix = [animal title_prefix_original];
        for chi = 1 : nchns -1
            chnnamei = chnsNames{chi};
            for chj = chi + 1 : nchns
                chnnamej = chnsNames{chj};
                sigciCoh_1pair = squeeze(ciCoh_tAOI(chi, chj, :, :));
                
                titlename = [title_prefix chnnamei '-' chnnamej '-' pdcond];
                
                % plot
                ifig = figure('Position', pos_ifig);
                plot_dynamic_ciCohHist_1pair(sigciCoh_1pair, f_selected, t_selected_AOI, time0name, ...
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
        clear title_prefix
        
        
        
        
        %%% plot sig dynamic ciCoh
        title_prefix = [animal title_prefix_sig];
        
        % extract sigCicoh_tAOI
        sigCicoh_tAOI = [];
        for ti = 1 : length(psedociCohs_AOI)
            psedociCohs_1t = psedociCohs_AOI{ti};
            ciCoh_tAOI_1t = squeeze(ciCoh_tAOI(:, :, :, ti));
            sigCicoh_1t = sigciCoh_extract(psedociCohs_1t, ciCoh_tAOI_1t, 'codesavefolder', savecodefolder);
            
            sigCicoh_tAOI = cat(4, sigCicoh_tAOI, sigCicoh_1t);
            clear tname psedociCohs_1t  ciCoh_tAOI_1t sigCicoh_1t
        end
        
        for chi = 1 : nchns -1
            chnnamei = chnsNames{chi};
            for chj = chi + 1 : nchns
                chnnamej = chnsNames{chj};
                sigciCoh_1pair = squeeze(sigCicoh_tAOI(chi, chj, :, :)); 
                
                titlename = [title_prefix chnnamei '-' chnnamej '-' pdcond];
                
                % plot
                ifig = figure('Position', pos_ifig);
                plot_dynamic_ciCohHist_1pair(sigciCoh_1pair, f_selected, t_selected_AOI, time0name, ...
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
        clear title_prefix
        
        
        clear pdcond 
        clear ciCohfile
        clear('ciCoh', 'psedociCohs', 'f_selected',  'T_chnsarea', 't_selected');
        clear ciCoh_tAOI  t_selected_AOI
    end
    clear animal input_folder 
end


function plot_dynamic_ciCohHist_1pair(ciCoh, f_selected, t_selected, time0name,...
    varargin)
%
%   Inputs:

% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'cbarTicks', [], @(x) assert(isvector(x) && isnumeric(x)));
addParameter(p, 'cbarStr', '', @isstr);
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_titlename', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'titlename', '', @isstr);

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig = p.Results.fig;
cbarTicks = p.Results.cbarTicks;
cbarStr = p.Results.cbarStr;
show_ylabel = p.Results.show_ylabel;
show_xlabel = p.Results.show_xlabel;
show_titlename = p.Results.show_titlename;
show_colorbar = p.Results.show_colorbar;
titlename = p.Results.titlename;


% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



%%%  plot 

if isempty(fig)
    width = p.Results.width;
    height = p.Results.height;
    fig = figure('Position', [50 50 width height]);
    set(fig, 'PaperUnits', 'points');
    clear width height
end

ax = axes(fig, 'Units', 'pixels');

imagesc(ax, t_selected, f_selected, ciCoh); hold on;
colormap(jet)
set(gca,'YDir','normal')
if ~isempty(histClim)
    set(gca,'CLim', histClim)
end


% set time0name
xtklabels = xticklabels();
xtklabels{cellfun(@(x) strcmp(x, '0'), xtklabels)} = time0name;
xticklabels(xtklabels);

% time0name line
plot([0 0], ylim, '--r', 'LineWidth',2)


c = colorbar;

if ~isempty(cbarStr)
    c.Label.String = cbarStr;
end
if ~isempty(cbarTicks)
    set(c, 'Ticks', cbarTicks, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roma')
end
c.Visible = 'off';


%%% show inf
if show_ylabel
    ylabel('Frequence (Hz)','FontName','Times New Roma')
end

if show_xlabel
    xlabel('time/s','FontName','Times New Roma')
end

if show_titlename && ~isempty(titlename)
    title(titlename, 'FontSize', 10, 'FontWeight', 'normal')
end

if show_colorbar
    c.Visible = 'on';
end



