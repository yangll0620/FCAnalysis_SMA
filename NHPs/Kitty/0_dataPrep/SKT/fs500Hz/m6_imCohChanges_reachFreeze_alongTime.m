function m6_imCohChanges_reachFreeze_alongTime(varargin)
% plot cicoh changes of early, middle, late and after200ms reach freeze
% relative to previous freeze phase
%
%   Example Usage:
%           m6_imCohChanges_reachFreeze_alongTime();
%           m6_imCohChanges_reachFreeze_alongTime('shuffleN_psedoTest', 500, 'newRun', true)
%           m6_imCohChanges_reachFreeze_alongTime('plotCiCohChanges', false, 'plotCiCoh', true)
%   
%   Input:
%       Name-Value:
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500
%           plotCiCohChanges - plot ciCohChanges true (default) or false
%           plotCiCoh - plot ciCoh true (default) or false

codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));


% parse params
p = inputParser;
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plotCiCohChanges', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plotCiCoh', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
plotCiCohChanges = p.Results.plotCiCohChanges;
plotCiCoh = p.Results.plotCiCoh;



% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


% animal
animal = animal_extract(codecorresfolder);


%%  input setup
lfpfile = fullfile(codecorresParentfolder, 'm5_Kitty_imCoh_reachFreeze', 'ciCohs-reachFreeze.mat');



%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code start here

load(lfpfile, 'lfpsegs_Freeze');
reachfreezeTypes = fieldnames(lfpsegs_Freeze.ReachFreeze);
clear lfpsegs_Freeze

nfreezeTypes = length(reachfreezeTypes);
ciCoh_Changes_file = fullfile(savefolder, 'ciCohsChanges-reachFreeze-basedPrePhases.mat');
if (~exist(ciCoh_Changes_file, 'file') || newRun)
    
    % load and save
    load(lfpfile, 'ciCohs', 'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'T_chnsarea', 'f_AOI','f_selected');
    save(ciCoh_Changes_file, 'ciCohs', 'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'f_AOI','T_chnsarea', 'f_selected');
    clear('ciCohs', 'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'T_chnsarea', 'f_AOI', 'f_selected');
    
    % calculate ciCohChanges
    load(lfpfile, 'ciCohs');
    for fri_base = 1 : nfreezeTypes - 1
        subfreezeType_base = reachfreezeTypes{fri_base};
        ciCoh_base = ciCohs.ReachFreeze.(subfreezeType_base);

        for fri_comp = fri_base + 1 : nfreezeTypes
            subfreezeType_comp = reachfreezeTypes{fri_comp};
            
            ciCoh_comp = ciCohs.ReachFreeze.(subfreezeType_comp);

            % calculate ciCohChangs
            ciCohChanges.ReachFreeze.(['b_' subfreezeType_base]).(subfreezeType_comp) = ciCoh_comp - ciCoh_base;

            clear subfreezeType ciCoh_comp
        end
        clear subfreezeType_base ciCoh_base
    end
    save(ciCoh_Changes_file, 'ciCohChanges', '-append');
end

% psedoCicohChanges
freezeType = 'ReachFreeze';
load(ciCoh_Changes_file, 'lfpsegs_Freeze', 'fs', 'f_AOI')
for fri_base = 1 : nfreezeTypes - 1
    subfreezeType_base = reachfreezeTypes{fri_base};
    lfpsegs_base = lfpsegs_Freeze.ReachFreeze.(subfreezeType_base);

    for fri_comp = 1 : nfreezeTypes
        subfreezeType_comp = reachfreezeTypes{fri_comp};
        lfpsegs_comp = lfpsegs_Freeze.ReachFreeze.(subfreezeType_comp);

        % psedo ciCohChanges
        psedoCiCohChanges_FreezeBasePrePhases_extract_save(shuffleN_psedoTest, lfpsegs_comp, lfpsegs_base, fs, f_AOI, ciCoh_Changes_file, ['b' subfreezeType_base], freezeType, subfreezeType_comp);

        clear subfreezeType lfpsegs_comp 
    end
    clear subfreezeType_base lfpsegs_base fri_comp
end
clear('lfpsegs_Freeze', 'fs', 'f_AOI')


% plot ciCohChanges
if plotCiCohChanges 
    load(ciCoh_Changes_file, 'ciCohChanges','psedociCohChanges', 'f_selected',  'T_chnsarea')
    
    cbarStr = 'ciCohChange';
    cbarTicks = [-1 0 1];
    histClim = [-1 1];
    
    for fri_base = 1 : nfreezeTypes - 1
        subfreezeType_base = reachfreezeTypes{fri_base};
    
        show_xticklabels = false;
        show_xlabel = false;
        if fri_base == nfreezeTypes - 1
            show_xticklabels = true;
            show_xlabel = true;
        end
        for fri_comp = fri_base + 1 : nfreezeTypes
            subfreezeType_comp = reachfreezeTypes{fri_comp};
    
            ciCohchange = ciCohChanges.ReachFreeze.(['b_' subfreezeType_base]).(subfreezeType_comp);
            psedoiCohChange = psedociCohChanges.ReachFreeze.(['b' subfreezeType_base]).(subfreezeType_comp);
    
            [sigciCohChanges]= sigciCoh_extract(psedoiCohChange, ciCohchange);
            [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
    
    
            show_titlename = false;
            show_yticklabels = false;
            show_colorbar = false;
            if fri_comp == 1
                show_yticklabels = true;
            end
            if fri_comp == length(reachfreezeTypes)
                show_colorbar = true;
            end
    
    
            % plot and save ciCoh Histogram image
            ifig = figure('Position', [50 50 400 200]);
            set(ifig, 'PaperUnits', 'points');
            plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, ['reachFreeze ' subfreezeType_comp '-' subfreezeType_base], 'histClim', histClim,...
                'codesavefolder', savecodefolder, 'cbarStr', cbarStr, 'cbarTicks', cbarTicks, ...
                'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
                'fig', ifig);
            subfilename = ['reachFreeze-base' subfreezeType_base '-' subfreezeType_comp]; 
            print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
            print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
            close(ifig)


            clear subfreezeType 
            clear ciCohchange psedoiCohChange sigciCohChanges sigciCohChanges_flatten chnPairNames
            clear show_titlename show_yticklabels show_colorbar
            clear ifig subfilename
    
        end
    
        clear eBasePhase show_xlabel show_xticklabels fri
    end
    clear('ciCohChanges','psedociCohChanges', 'f_selected',  'T_chnsarea')
end


% plot ciCoh
if plotCiCoh
    load(ciCoh_Changes_file, 'ciCohs','psedociCohs', 'f_selected',  'T_chnsarea');

    cbarStr = 'ciCoh';
    cbarTicks = [0 0.5 1];
    histClim = [0 1];

    show_xticklabels = false;
    show_xlabel = false;

    for fri = 1 : length(reachfreezeTypes)
        subfreezeType_comp = reachfreezeTypes{fri};

        ciCoh = ciCohs.ReachFreeze.(subfreezeType_comp);
        psedociCoh = psedociCohs.ReachFreeze.(subfreezeType_comp);

        [sigciCoh]= sigciCoh_extract(psedociCoh, ciCoh);
        [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);


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
        ifig = figure('Position', [50 50 400 200]);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram2(sigciCoh_flatten, chnPairNames, f_selected, ['reachFreeze ' subfreezeType_comp], ...
            'histClim', histClim, 'cbarTicks', cbarTicks, 'cbarStr', cbarStr, ...
            'codesavefolder', savecodefolder, ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig);
        subfilename = ['reachFreeze-' subfreezeType_comp];
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
        close(ifig)

        clear subfreezeType 
        clear ciCoh psedociCoh sigciCoh sigciCoh_flatten chnPairNames
        clear show_titlename show_yticklabels show_colorbar
        clear ifig subfilename
    end

    clear show_xlabel show_xticklabels fri
    clear('ciCohs','psedociCohs', 'f_selected',  'T_chnsarea')
end


function psedoCiCohChanges_FreezeBasePrePhases_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile,  eBasePhase, freezeType, subfreezeType)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohPhasefile:
%       
%   
% save psedoiCohChanges: nchns * nchns * nf * nshuffle, saved to ciCohChangesfile



load(ciCohChangesfile, 'psedociCohChanges');

if(~exist('psedociCohChanges', 'var'))
    psedociCohChanges = struct();
end

if ~isfield(psedociCohChanges, freezeType) || ~isfield(psedociCohChanges.(freezeType), eBasePhase) || ~isfield(psedociCohChanges.(freezeType).(eBasePhase), subfreezeType)
    psedociCohChanges.(freezeType).(eBasePhase).(subfreezeType) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges.(freezeType).(eBasePhase).(subfreezeType), 4) + 1;
end
lfp_combined = cat(3, lfptrials_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_comp = lfp_combined(:, :, masksBase);
    psedolfp_Base = lfp_combined(:, :, ~masksBase);
    
    [~, psedoiCoh_comp, ~] = ciCoh_trialDeltaPhi(psedolfp_comp, fs, f_AOI);
    [~, psedoiCoh_Base, ~] = ciCoh_trialDeltaPhi(psedolfp_Base, fs, f_AOI);
    
    psedociCohChanges.(freezeType).(eBasePhase).(subfreezeType) = cat(4, psedociCohChanges.(freezeType).(eBasePhase).(subfreezeType), psedoiCoh_comp - psedoiCoh_Base);
    
    if(mod(si, 100) == 0)
        disp([freezeType '-' eBasePhase '-' subfreezeType ' pesdo ciCoh Changes test ' num2str(si)])
        save(ciCohChangesfile, 'psedociCohChanges', '-append');
    end
    
    clear randomSKTInds randomRestInds psedolfp_comp psedoiCoh_rest
    clear psedoiCoh_comp psedoiCoh_rest
end

function plot_ciCohHistogram2(ciCoh_flatten, chnPairNames, f_selected, titlename, varargin)
%
%   Inputs:
%       ciCoh_flatten:
%       chnPairNames
%       f_selected
%       titlename
%       histClim
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
%           'histClim' - ciCoh histogram clim (default [0 1])
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'cbarTicks' - vector, default [0 0.5 1]
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_titlename' - show (true, default) or not show (false) titlename
%           'show_colorbar' - show (true, default) or not show (false) colorbar
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default []
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default []
%           'fontsize1' - font size for title, default 12
%           'fontsize2' - font size for , default 10
%           'fontname' - font name , default 'Times New Roman'



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [0 1], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'cbarTicks', [0 0.5 1], @(x) assert(isvector(x) && isnumeric(x)));
addParameter(p, 'cbarStr', 'ciCoh', @isstr);
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_titlename', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'innerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'fontsize1', 11, @(x) assert(isscalar(x) && isnumeric(x)));
addParameter(p, 'fontsize2', 10, @(x) assert(isscalar(x) && isnumeric(x)));
addParameter(p, 'fontname', 'Times New Roman', @isstr);
addParameter(p, 'width', 560, @(x) assert(isscalar(x) && isnumeric(x)));
addParameter(p, 'height', 420, @(x) assert(isscalar(x) && isnumeric(x)));


parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig = p.Results.fig;
cbarTicks = p.Results.cbarTicks;
cbarStr = p.Results.cbarStr;
show_xticklabels = p.Results.show_xticklabels;
show_yticklabels = p.Results.show_yticklabels;
show_xlabel = p.Results.show_xlabel;
show_titlename = p.Results.show_titlename;
show_colorbar = p.Results.show_colorbar;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin; 
fontsize1 = p.Results.fontsize1;
fontsize2 = p.Results.fontsize2;
fontname = p.Results.fontname;



% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% plot
if isempty(fig)
    width = p.Results.width;
    height = p.Results.height;
    fig = figure('Position', [50 50 width height]);
    set(fig, 'PaperUnits', 'points');
    clear width height
end
ax = axes(fig, 'Units', 'pixels');
imagesc(ax, ciCoh_flatten)
colormap(jet)
set(gca,'CLim', histClim)

% set axes position
if ~isempty(outerposMargin) && ~isempty(innerposMargin)
    fig_pos = fig.Position;
    fig_width = fig_pos(3);
    fig_height = fig_pos(4);
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
end

%%% plot the line to separete the pairs
chnPair_prev = '';
for ci = 1: length(chnPairNames)
    chnPair = chnPairNames{ci};
    
    % replace M1-stn0-1 to M1-STN
    s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
    if ~isempty(s_stn)
        for si = 1 : length(s_stn)
            chnPair = strrep(chnPair, s_stn{si}, 'STN');
        end
    end
    % replace M1-stn0-1 to M1-STN
    s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
    if ~isempty(s_gp)
        for si = 1 : length(s_gp)
            chnPair = strrep(chnPair, s_gp{si}, 'GP');
        end
    end
    
    if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
        hold on; plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'w--')
        % Create line
    end
    chnPair_prev = chnPair;
    chnPairNames{ci} = chnPair;
    
    
    clear s_stn s_gp chnPair
end

c = colorbar;
c.Label.String = cbarStr;
if ~isempty(cbarTicks)
    set(c, 'Ticks', cbarTicks, 'FontSize', fontsize1, 'FontWeight', 'bold', 'FontName', fontname)
end
c.Visible = 'off';


%%% show inf
[npairs, nf] = size(ciCoh_flatten);
if show_xticklabels
    xticks([1:nf])
    xticklabels(round(f_selected))  
    set(gca, 'fontsize',fontsize2, 'FontName', fontname, 'FontWeight', 'bold')
else
    xticks([]);
end

if show_yticklabels
    yticks([1:npairs]);
    set(gca,'YTickLabel',chnPairNames,'fontsize',11, 'FontName', fontname, 'FontWeight', 'bold')
else 
    yticks([]);
end

if show_xlabel
    xlabel('Frequency (Hz)', 'fontsize', 12, 'FontName', fontname, 'FontWeight', 'bold')
end
if show_titlename
    title(titlename, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_colorbar
    c.Visible = 'on';
end
