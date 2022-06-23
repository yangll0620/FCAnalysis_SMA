function m4_imCohChanges_phases2preMove_perslowFast(varargin)
% cicoh changes for early and late reach separete slow and fast reach
% cases
%
%   Example Usage:
%           m4_imCohChanges_phases2preMove_perslowFast('t_slowFastTheshold', 1.5, 'shuffleN_psedoTest', 500, 'newRun', true)
%   
%   Input:
%
%       Name-Value:
%           t_slowFastTheshold - time theshold for slow and fast reach
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
addParameter(p, 't_slowFastTheshold', 1.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
t_slowFastTheshold = p.Results.t_slowFastTheshold;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


% animal
animal = animal_extract(codecorresfolder);


%%  input setup
inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');



%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


f_AOI = [8 40];



%% Code start here
pdcond = 'moderate';
eCompPhases = {'earlyReach';  'PeakV'; 'lateReach'};
eBasePhase = 'preMove';

lfpfiles = dir(fullfile(inputfolder, ['*' pdcond '_*.mat']));

ciCoh_Changes_file = fullfile(savefolder, 'ciCohs-Changes.mat');


if (~exist(ciCoh_Changes_file, 'file') || newRun)

    [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI_extract(eBasePhase, animal, pdcond, 'codesavefolder', savecodefolder);

    % extract lfptrials of slow and fast trials for base phase
    [baselfptrials_fast, fs, T_chnsarea] = lfptrials_animalK_align2(lfpfiles, align2, t_AOI,...
        't_min_reach', 0.5, 't_max_reach', t_slowFastTheshold);
    [baselfptrials_slow, ~, ~] = lfptrials_animalK_align2(lfpfiles, align2, t_AOI, ...
        't_min_reach', t_slowFastTheshold, 't_max_reach', inf);

    lfptrials.fastReach.(['b_' eBasePhase]) = baselfptrials_fast;
    lfptrials.slowReach.(['b_' eBasePhase]) = baselfptrials_slow;


    % calculate base ciCohs for fast and slow cases
    [~, baseciCoh_fast, f_selected]= ciCoh_trialDeltaPhi(baselfptrials_fast, fs, f_AOI);
    [~, baseciCoh_slow, ~]= ciCoh_trialDeltaPhi(baselfptrials_slow, fs, f_AOI);
    ciCohs.fastReach.(['b_' eBasePhase]) = baseciCoh_fast;
    ciCohs.slowReach.(['b_' eBasePhase]) = baseciCoh_slow;


    for ei = 1 : length(eCompPhases)
        eCompPhase = eCompPhases{ei};
        [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI_extract(eCompPhase, animal, pdcond, 'codesavefolder', savecodefolder);

        % extract lfptrials of slow and fast trials
        [complfptrials_fast, ~, ~] = lfptrials_animalK_align2(lfpfiles, align2, t_AOI,...
            't_min_reach', 0.5, 't_max_reach', t_slowFastTheshold);
        [complfptrials_slow, ~, ~] = lfptrials_animalK_align2(lfpfiles, align2, t_AOI, ...
            't_min_reach', t_slowFastTheshold, 't_max_reach', inf);


        lfptrials.fastReach.(eCompPhase) = complfptrials_fast;
        lfptrials.slowReach.(eCompPhase) = complfptrials_slow;


        % calculate compare ciCohs for fast and slow cases
        [~, compciCoh_fast, ~]= ciCoh_trialDeltaPhi(complfptrials_fast, fs, f_AOI);
        [~, compciCoh_slow, ~]= ciCoh_trialDeltaPhi(complfptrials_slow, fs, f_AOI);

        ciCohs.fastReach.(eCompPhase) = compciCoh_fast;
        ciCohs.slowReach.(eCompPhase) = compciCoh_slow;


        % calculate ciCohChangs
        ciCohChanges.fastReach.(['b_' eBasePhase]).(eCompPhase) = compciCoh_fast - baseciCoh_fast;
        ciCohChanges.slowReach.(['b_' eBasePhase]).(eCompPhase) = compciCoh_slow - baseciCoh_slow;
        
        clear eCompPhase align2 t_AOI
        clear complfptrials_fast complfptrials_slow
        clear compciCoh_fast compciCoh_slow
    end

    save(ciCoh_Changes_file, 'lfptrials', 'fs', 'T_chnsarea', 'f_AOI');
    save(ciCoh_Changes_file, 'ciCohs', 'f_selected', 'ciCohChanges', '-append');

    clear align2 t_AOI fs T_chnsarea 
    clear baselfptrials_fast baselfptrials_slow baseciCoh_fast baseciCoh_slow 
    clear('lfptrials', 'fs', 'T_chnsarea', 'f_AOI', 'ciCohs', 'ciCohChanges', 'f_selected')
end


% psedoCicoh
load(ciCoh_Changes_file, 'lfptrials', 'fs', 'f_AOI')
for ei = 1 : length(eCompPhases)
    eCompPhase = eCompPhases{ei};


    baselfptrials_fast = lfptrials.fastReach.(['b_' eBasePhase]);
    baselfptrials_slow = lfptrials.slowReach.(['b_' eBasePhase]);

    complfptrials_fast = lfptrials.fastReach.(eCompPhase);
    complfptrials_slow = lfptrials.slowReach.(eCompPhase);

    psedociCohChanges_ePhaseComp_extract_save(shuffleN_psedoTest, complfptrials_fast, baselfptrials_fast, fs, f_AOI, ciCoh_Changes_file,  'fastReach', ['b_' eBasePhase], eCompPhase);
    psedociCohChanges_ePhaseComp_extract_save(shuffleN_psedoTest, complfptrials_slow, baselfptrials_slow, fs, f_AOI, ciCoh_Changes_file,  'slowReach', ['b_' eBasePhase], eCompPhase);

    clear eCompPhase lfptrials_slow lfptrials_fast
end
clear('lfptrials', 'fs', 'f_AOI')


% plot
slowFastTypes = {'fastReach';'slowReach'};
load(ciCoh_Changes_file, 'ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea')

for sfi = 1 : length(slowFastTypes)
    slowFastType = slowFastTypes{sfi};

    show_yticklabels = false;
    show_colorbar = false;
    if sfi == 1
        show_yticklabels = true;
    end
    if sfi == length(slowFastTypes)
        show_colorbar = true;
    end

    for ei = 1 : length(eCompPhases)
        eCompPhase = eCompPhases{ei};

        ciCohchange = ciCohChanges.(slowFastType).(['b_' eBasePhase]).(eCompPhase);
        psedoiCohChange = psedociCohChanges.(slowFastType).(['b_' eBasePhase]).(eCompPhase);

        [sigciCohChanges]= sigciCoh_extract(psedoiCohChange, ciCohchange);
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);


        show_xticklabels = false;
        show_xlabel = false;
        show_titlename = false;
        if ei == length(eCompPhases)
            show_xticklabels = true;
            show_xlabel = true;
        end

        % plot and save ciCoh Histogram image
        ifig = figure('Position', [50 50 400 200]);
        set(ifig, 'PaperUnits', 'points');
        plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, [slowFastType '-b_' eBasePhase '-'  eCompPhase], 'histClim', [-1 1],...
            'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
            'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
            'fig', ifig);

        subfilename = [slowFastType 'Trials-b' eBasePhase '-'  eCompPhase]; % 'Jo-mild-Bnormal-preMove'
        print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
        print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
        close(ifig)

        clear eCompPhase ciCohchange

    end

    clear slowFastType
end


function psedociCohChanges_ePhaseComp_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile,  slowFastType, eBasePhase, eCompPhase)
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

if ~isfield(psedociCohChanges, slowFastType) || ~isfield(psedociCohChanges.(slowFastType), eBasePhase) || ~isfield(psedociCohChanges.(slowFastType).(eBasePhase), eCompPhase)
    psedociCohChanges.(slowFastType).(eBasePhase).(eCompPhase) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges.(slowFastType).(eBasePhase).(eCompPhase), 4) + 1;
end

lfp_combined = cat(3, lfptrials_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
  
    psedolfp_Base = lfp_combined(:, :, masksBase);
    psedolfp_comp = lfp_combined(:, :, ~masksBase);

    [~, psedoiCoh_Base, ~] = ciCoh_trialDeltaPhi(psedolfp_Base, fs, f_AOI);
    [~, psedoiCoh_comp, ~] = ciCoh_trialDeltaPhi(psedolfp_comp, fs, f_AOI);
    
    psedociCohChanges.(slowFastType).(eBasePhase).(eCompPhase) = cat(4, psedociCohChanges.(slowFastType).(eBasePhase).(eCompPhase), psedoiCoh_comp - psedoiCoh_Base);
    
    if(mod(si, 100) == 0)
        save(ciCohChangesfile, 'psedociCohChanges', '-append');
        disp([slowFastType '-' eBasePhase '-' eCompPhase ' pesdociCohChanges test ' num2str(si)])
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




function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_animalK_align2(files, align2, tdur_trial, varargin)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately
% [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2PeakV(files, [t_AOI(1) t_AOI(2)], 'codesavefolder', savecodefolder);
%
%   not include trials with t_reach <0.2s
% 
%         Args:
%             align2: the event to be aligned 
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])     
%
%             
%
%       Name-Value: 
%            't_min_reach' - the minimal reach time used, default 0.5 s
%
%   
%           't_max_reach' - the max reach time used, default inf
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


% parse params
p = inputParser;
addParameter(p, 't_min_reach', 0.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 't_max_reach', inf, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
t_min_reach = p.Results.t_min_reach;
t_max_reach = p.Results.t_max_reach;


% code start here
coli_align2 = uint32(align2);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);


load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'selectedTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    ntrials = length(lfpdata);
    for tri = 1: ntrials
        
        % ignore trials marked with 0
        if ~selectedTrials(tri)
            continue
        end
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_min_reach  || t_reach > t_max_reach
            clear t_reach
            continue
        end
        
        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
            idx_peakV_ma = idx + idx_reachonset_ma -1;
            t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
            
            if t_reachonset2peakV < 0.3 || t_peakV2reach < 0.3
                clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
                clear t_reachonset2peakV  t_peakV2reach
                continue;
            end
            
            % extract trial with t_dur
            idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
            idx_time0 = idx_peakV_lfp;
            
            clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
            clear t_reachonset2peakV  t_peakV2reach
            clear idx_peakV_lfp
        else
            idx_time0 = T_idxevent_lfp{tri, coli_align2}; 
        end
        
        
        % extract phase for 1 trial
        lfp_1trial = lfpdata{tri};
        idxdur = round(tdur_trial * fs_lfp) + idx_time0;
        if idxdur(1) == 0
            idxdur(1) = 1;
        else
            idxdur(1) = idxdur(1) + 1;
        end
        lfp_phase_1trial = lfp_1trial(:,idxdur(1) :idxdur(2));
           
        % cat into lfptrials
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach idxdur lfp_phase_1trial lfp_1trial
    end
end


