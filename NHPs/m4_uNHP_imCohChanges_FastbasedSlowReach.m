function m4_uNHP_imCohChanges_FastbasedSlowReach(animal, varargin)
% plot cicoh changes of slow trials compared to the fast trials
%
%   Example Usage:
%           m4_uNHP_imCohChanges_FastbasedSlowReach('Jo', 'nslowfast_J', 30, 'shuffleN_psedoTest', 500, 'newRun', true)
%   
%   Input:
%
%       Name-Value:
%           nslowfast_J - fast and slow ntrials used for animal J default 30
%           nslowfast_K - fast and slow ntrials used for animal K, default 30
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
addParameter(p, 'nslowfast_J', 30, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'nslowfast_K', 30, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
nslowfast_J = p.Results.nslowfast_J;
nslowfast_K = p.Results.nslowfast_K;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;



%%  input setup
if strcmpi(animal, 'Kitty')
    lfpfolder = 'm2_segSKTData_SelectTrials_chnOfI';
    nslowfast = nslowfast_K;
end
if strcmpi(animal, 'Jo')
    lfpfolder = 'm2_SKTData_SelectTrials';
    nslowfast = nslowfast_J;
end

[~, ~, pipelinefolder, ~] = exp_subfolders();
inputfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', lfpfolder);


%% save setup
[~,codefilename,~] = fileparts(codefilepath); 
savecodefolder = fullfile(codefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', codefilename);
[savefolder, ~] = code_corresfolder(savecodefolder, false, false);

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);



%% Code start here
pdcond = 'mild';
f_AOI = [8 40];

ePhases = {'preMove'; 'earlyReach';  'PeakV'; 'lateReach'};

lfpfiles = dir(fullfile(inputfolder, ['*' pdcond '_*.mat']));
t_reachs = extract_ReachTime(lfpfiles);
t_reachs = sort(t_reachs);

t_slowReach_min = t_reachs(end-nslowfast + 1);
t_slowReach_max = t_reachs(end);
t_fastReach_min = t_reachs(1);
t_fastReach_max = t_reachs(nslowfast);

ciCoh_Changes_file = fullfile(savefolder, [animal 'ciCohChanges_FastbasedSlowReach_' pdcond '_ntrials'  num2str(nslowfast) '.mat']);


if (~exist(ciCoh_Changes_file, 'file') || newRun)

    for ei = 1 : length(ePhases)
        ePhase = ePhases{ei};
        [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI_extract(ePhase, animal, pdcond, 'codesavefolder', savecodefolder);

        % extract lfptrials of slow and fast trials
        if strcmpi(animal, 'Kitty')
            [lfptrials_fast, fs, T_chnsarea] = lfptrials_animalK_align2(lfpfiles, align2, t_AOI,...
                't_min_reach', t_fastReach_min, 't_max_reach', t_fastReach_max);
            [lfptrials_slow, ~, ~] = lfptrials_animalK_align2(lfpfiles, align2, t_AOI, ...
                't_min_reach', t_slowReach_min, 't_max_reach', t_slowReach_max);
        end

        if strcmpi(animal, 'Jo')
            [lfptrials_fast, fs, T_chnsarea] = lfptrials_animalJ_align2(lfpfiles, align2, t_AOI,...
                't_min_reach', t_fastReach_min, 't_max_reach', t_fastReach_max);
            [lfptrials_slow, ~, ~] = lfptrials_animalJ_align2(lfpfiles, align2, t_AOI, ...
                't_min_reach', t_slowReach_min, 't_max_reach', t_slowReach_max);
            
            chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);
            lfptrials_fast = lfptrials_fast(mask_chnOfI, :, :);
            lfptrials_slow = lfptrials_slow(mask_chnOfI, :, :);
            clear mask_chnOfI chnsOfI
        end


        lfptrials.slowReach.(ePhase) = lfptrials_slow;
        lfptrials.fastReach.(ePhase) = lfptrials_fast;


        % calculate ciCohs
        [~, ciCoh_slow, f_selected]= ciCoh_trialDeltaPhi(lfptrials_slow, fs, f_AOI);
        [~, ciCoh_fast, ~]= ciCoh_trialDeltaPhi(lfptrials_fast, fs, f_AOI);
        ciCohs.slowReach.(ePhase) = ciCoh_slow;
        ciCohs.fastReach.(ePhase) = ciCoh_fast;


        % calculate ciCohChangs
        ciCohChanges.Fast2Slow.(ePhase) = ciCoh_fast - ciCoh_slow;
        
        clear lfptrials_slow lfptrials_fast
        clear ciCoh_slow ciCoh_fast
    end

    save(ciCoh_Changes_file, 'lfptrials', 'fs', 'T_chnsarea', 'f_AOI');
    save(ciCoh_Changes_file, 'ciCohs', 'f_selected', 'ciCohChanges', '-append');

    clear fs T_chnsarea 
    clear align2name align2 t_AOI

end
clear('lfptrials', 'fs', 'T_chnsarea', 'f_AOI', 'ciCohs', 'ciCohChanges', 'f_selected')


% psedoCicoh
load(ciCoh_Changes_file, 'lfptrials', 'fs', 'f_AOI')
for ei = 1 : length(ePhases)
    ePhase = ePhases{ei};

    lfptrials_slow = lfptrials.slowReach.(ePhase);
    lfptrials_fast = lfptrials.fastReach.(ePhase);

    % psedo ciCoh for slow and fast individually
    psedoSlowFastCiCoh_extract_save(shuffleN_psedoTest, lfptrials_slow, fs, f_AOI, ciCoh_Changes_file, 'slowReach', ePhase);
    psedoSlowFastCiCoh_extract_save(shuffleN_psedoTest, lfptrials_fast, fs, f_AOI, ciCoh_Changes_file, 'fastReach', ePhase);

    % psedo ciCohChanges
    FastSlowBaseTye = 'Fast2Slow';
    psedoSlowFastciCohChanges_extract_save(shuffleN_psedoTest, lfptrials_fast, lfptrials_slow, fs, f_AOI, ciCoh_Changes_file, FastSlowBaseTye, ePhase);

    clear ePhase lfptrials_slow lfptrials_fast
end
clear('lfptrials', 'fs', 'f_AOI')


% plot
load(ciCoh_Changes_file, 'ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea')
for ei = 1 : length(ePhases)
    ePhase = ePhases{ei};

    ciCohchange = ciCohChanges.(ePhase);
    psedoiCohChange = psedociCohChanges.(ePhase);

    [sigciCohChanges]= sigciCoh_extract(psedoiCohChange, ciCohchange);
    [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);


    show_yticklabels = false;
    show_colorbar = true;
    show_xlabel = false;
    show_xticklabels = false;
    show_titlename = false;
    if ei == length(ePhases)
        show_xlabel = true;
        show_xticklabels = true;
    end


    % plot and save ciCoh Histogram image
    ifig = figure('Position', [50 50 400 200]);
    set(ifig, 'PaperUnits', 'points');
    plot_ciCohHistogram2(sigciCohChanges_flatten, chnPairNames, f_selected, [ePhase ' slowReach-fastReach'], 'histClim', [-1 1],...
        'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1], ...
        'show_xticklabels', show_xticklabels, 'show_yticklabels', show_yticklabels, 'show_xlabel', show_xlabel, 'show_titlename', show_titlename,'show_colorbar', show_colorbar, ...
        'fig', ifig);
    subfilename = ['ciCohChanges-' 'FastbasedSlowReach-' ePhase]; % 'Jo-mild-Bnormal-preMove'
    print(ifig, fullfile(savefolder, subfilename), '-painters', '-depsc')
    print(ifig, fullfile(savefolder, subfilename), '-dpng', '-r1000')
    close(ifig)

    clear ePhase ciCohchange
end


function psedoSlowFastciCohChanges_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile, FastSlowBaseTye, ePhase)
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

if (~isfield(psedociCohChanges, FastSlowBaseTye)) || (~isfield(psedociCohChanges.(FastSlowBaseTye), ePhase))
    psedociCohChanges.(FastSlowBaseTye).(ePhase) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges.(FastSlowBaseTye).(ePhase), 4) + 1;
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
    
    [~, psedoiCoh_comp, ~] = ciCoh_trialDeltaPhi(psedolfp_comp, fs, f_AOI);
    [~, psedoiCoh_Base, ~] = ciCoh_trialDeltaPhi(psedolfp_Base, fs, f_AOI);
    
    psedociCohChanges.(FastSlowBaseTye).(ePhase) = cat(4, psedociCohChanges.(FastSlowBaseTye).(ePhase), psedoiCoh_comp - psedoiCoh_Base);
    
    if(mod(si, 100) == 0)
        disp([FastSlowBaseTye '-' ePhase ' pesdo ciCoh Changes test ' num2str(si)])
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




function psedoSlowFastCiCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohfile, slowFastType, ePhase)
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

if ~isfield(psedociCohs, slowFastType) || ~isfield(psedociCohs.(slowFastType), ePhase) 
    psedociCohs.(slowFastType).(ePhase) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs.(slowFastType).(ePhase), 4) + 1;
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
    psedociCohs.(slowFastType).(ePhase) = cat(4, psedociCohs.(slowFastType).(ePhase), psedociCohs1Time);
    clear psedociCohs1Time
    if(mod(si,100) == 0)
        disp([slowFastType '- ' ePhase '- psedo test ' num2str(si) ' times'])
        save(ciCohfile, 'psedociCohs', '-append');
    end
end


function t_reachs = extract_ReachTime(files)

coli_reachonset = 2;
coli_touch = 3;

t_event = [];
for fi = 1: length(files)
    file = fullfile(files(fi).folder, files(fi).name);
    load(file, 'T_idxevent_lfp', 'fs_lfp');
    var = whos('-file',file, 'goodTrials', 'selectedTrials');
    load(file, var.name);
    if strcmp(var.name, 'selectedTrials')
        goodTrials = selectedTrials;
        clear selectedTrials
    end
    t_event = cat(1, t_event, T_idxevent_lfp{goodTrials, :}/ fs_lfp);
    clear T_idxevent_lfp fs_lfp goodTrials
end
t_reachs= t_event(:, coli_touch) - t_event(:, coli_reachonset);
clear pdcond files t_event fi



function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_animalK_align2(files, align2, tdur_trial, varargin)
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

function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_animalJ_align2(files, align2, tdur_trial, varargin)
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
    file = fullfile(files(filei).folder, filename);
    load(file, 'lfpdata', 'T_idxevent_lfp', 'goodTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
    selectedTrials = goodTrials;
    clear file goodTrials
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    ntrials = length(selectedTrials);
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
            [~, idx] = max(smoothWspeed_trial(idx_reachonset_ma: idx_reach_ma, tri));
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
        lfp_1trial = lfpdata(:, :, tri);
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

