function m5_imCohChanges_Freeze2ReachPhases(varargin)
% plot cicoh changes of early, middle, late and after200ms all freeze (combine Init, Reach and Mani Freezes)
% relative to premove and early reach
%
%   Example Usage:
%           m5_imCohChanges_Freeze2ReachPhases('shuffleN_psedoTest', 500, 'newRun', true)
%   
%   Input:
%       Name-Value:
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500

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

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


% animal
animal = animal_extract(codecorresfolder);


%%  input setup

lfpfile = fullfile(codecorresParentfolder, 'm4_fs500Hz_FreezeSegs_extract', 'Kitty-FreezeSegs-ReachTrials-FreezeSegs.mat');

f_AOI = [8 40];


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code start here

load(lfpfile, 'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'T_chnsarea');

eBasePhases = {'preMove'; 'earlyReach'};
reachfreezeTypes = fieldnames(lfpsegs_Freeze.ReachFreeze);

ciCoh_Changes_file = fullfile(savefolder, 'ciCohs-Changes-Freeze2reachPhases.mat');
if (~exist(ciCoh_Changes_file, 'file') || newRun)

    for ebi = 1 : length(eBasePhases)
        eBasePhase = eBasePhases{ebi};

        lfpbase = lfptrials_Reach.(eBasePhase);
        [~, ciCoh_base, f_selected]= ciCoh_trialDeltaPhi(lfpbase, fs, f_AOI);

        ciCohs.(eBasePhase) = ciCoh_base;

        for fri = 1 : length(reachfreezeTypes)
            subfreezeType = reachfreezeTypes{fri};

            % combined all freezeTypes
            lfpsegs = cat(3, lfpsegs_Freeze.InitFreeze.(subfreezeType), lfpsegs_Freeze.ReachFreeze.(subfreezeType), lfpsegs_Freeze.ManiFreeze.(subfreezeType));

            [~, ciCoh_segs, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);

            ciCohs.(subfreezeType) = ciCoh_segs;


            % calculate ciCohChangs
            ciCohChanges.(['b_' eBasePhase]).(subfreezeType) = ciCoh_segs - ciCoh_base;


            clear subfreezeType lfpsegs ciCoh_segs
        end
        clear eBasePhase lfpbase ciCoh_base
    end
    save(ciCoh_Changes_file, 'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'f_AOI','T_chnsarea');
    save(ciCoh_Changes_file, 'ciCohs', 'f_selected', 'ciCohChanges', '-append');

    clear('lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'T_chnsarea');
    clear('ciCohs', 'f_selected', 'ciCohChanges');
end

% psedoCicoh and psedoCicohChanges
load(ciCoh_Changes_file, 'lfpsegs_Freeze', 'lfptrials_Reach', 'fs', 'f_AOI')
for ebi = 1 : length(eBasePhases)
    eBasePhase = eBasePhases{ebi};

    lfpbase = lfptrials_Reach.(eBasePhase);

    % base psedo Cicoh 
    psedoFreezeCiCoh_extract_save(shuffleN_psedoTest, lfpbase, fs, f_AOI, ciCoh_Changes_file, eBasePhase);

    for fri = 1: length(reachfreezeTypes)
        subfreezeType = reachfreezeTypes{fri};

        % combined all freezeTypes
        lfpsegs = cat(3, lfpsegs_Freeze.InitFreeze.(subfreezeType), lfpsegs_Freeze.ReachFreeze.(subfreezeType), lfpsegs_Freeze.ManiFreeze.(subfreezeType));

        % psedo Cicoh
        psedoFreezeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCoh_Changes_file, subfreezeType);


        % psedo ciCohChanges
        psedoFreezeReachCiCohChanges_extract_save(shuffleN_psedoTest, lfpsegs, lfpbase, fs, f_AOI, ciCoh_Changes_file, ['b' eBasePhase], subfreezeType);

        clear subfreezeType lfpsegs 
    end


    clear eBasePhase lfpbase fri
end
clear('lfpsegs_Freeze', 'lfptrials_Reach', 'fs', 'f_AOI')


% plot
load(ciCoh_Changes_file, 'ciCohChanges','psedociCohChanges', 'f_selected',  'T_chnsarea')
for ebi = 1 : length(eBasePhases)
    eBasePhase = eBasePhases{ebi};

    show_xticklabels = false;
    show_xlabel = false;
    if ebi == length(eBasePhases)
        show_xticklabels = true;
        show_xlabel = true;
    end
    for fri = 1 : length(reachfreezeTypes)
        subfreezeType = reachfreezeTypes{fri};


        ciCohchange = ciCohChanges.(['b_' eBasePhase]).(subfreezeType);
        psedoiCohChange = psedociCohChanges.(['b' eBasePhase]).(subfreezeType);

        [sigciCohChanges]= sigciCoh_extract(psedoiCohChange, ciCohchange);
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
        ifig = figure('Position', [50 50 400 200]);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, ['reachFreeze ' subfreezeType '-' eBasePhase], 'histClim', [-1 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig);
        subfilename = ['Freeze2ReachPhase-' eBasePhase '-' subfreezeType]; 
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
        close(ifig)

    end
end




function psedoFreezeCiCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohfile,  subType)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohfile:
%       ePhase: 
%       slowFastType: 'slowReach' or 'fastReach'
%             
%   
% save psedociCohs.(slowFastType).(ePhase): nchns * nchns * nf * nshuffle, saved to ciCohfile


nchns = size(lfptrials, 1);

load(ciCohfile, 'psedociCohs');
if ~exist('psedociCohs', 'var')
    psedociCohs = struct();
end

if ~isfield(psedociCohs, subType)
    psedociCohs.(subType) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs.(subType), 4) + 1;
end

for si = shuffi_str : suffi_end
    for chni = 1 : nchns -1
        lfp1 = squeeze(lfptrials(chni, :, :));
        for chnj = chni + 1 : nchns
            lfp2 = squeeze(lfptrials(chnj, :, :));
            [psedociCoh, ~] = psedo_ciCohSKT_FFT_NoAmp(lfp1, lfp2, fs, f_AOI);
            
            if chni == 1 && chnj == chni + 1
                nf = size(psedociCoh, 1);
                psedociCohs1Time = zeros(nchns, nchns, nf, 1);
                clear nf
            end
            
            psedociCohs1Time(chni, chnj, :, :) = psedociCoh;
            clear lfp2 psedociCoh
        end
        clear lfp1
    end
    psedociCohs.(subType) = cat(4, psedociCohs.(subType), psedociCohs1Time);
    clear psedociCohs1Time
    if(mod(si,100) == 0)
        save(ciCohfile, 'psedociCohs', '-append');
        disp([subType '- psedo test ' num2str(si) ' times'])
    end
end


function psedoFreezeReachCiCohChanges_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile,  eBasePhase, subfreezeType)
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

if ~isfield(psedociCohChanges, eBasePhase) || ~isfield(psedociCohChanges.(eBasePhase), subfreezeType) 
    psedociCohChanges.(eBasePhase).(subfreezeType) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges.(eBasePhase).(subfreezeType), 4) + 1;
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
    
    psedociCohChanges.(eBasePhase).(subfreezeType) = cat(4, psedociCohChanges.(eBasePhase).(subfreezeType), psedoiCoh_comp - psedoiCoh_Base);
    
    if si == shuffi_str
        disp([eBasePhase  '-' subfreezeType ' pesdociCohChanges Start ']) 
    end

    if(mod(si, 40) == 0)
        save(ciCohChangesfile, 'psedociCohChanges', '-append');
        disp([eBasePhase  '-' subfreezeType ' pesdo ciCoh Changes test ' num2str(si)])   
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
