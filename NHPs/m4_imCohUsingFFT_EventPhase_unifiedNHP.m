function m4_imCohUsingFFT_EventPhase_unifiedNHP(varargin)
if nargin < 1
    animal = 'Jo';
else
    animal = varargin{1};
end

codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));

% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
segVFolder = false;
if strcmpi(animal, 'Kitty')
    segVFolder = true;
end
    
if segVFolder  
    SKTSubfolder = 'SKT_SegV';
    SKTDataSubfolder = 'm2_segSKTData_SelectTrials_goodReach';
else
    SKTSubfolder = 'SKT';
    SKTDataSubfolder = 'm2_SKTData_SelectTrials';
end
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);


%% save setup
savefolder = codecorresfolder;
copyfile2folder(codefilepath, savefolder);

ciCohPhasefile_prefix =[animal 'ciCohPhasefile'];

%%  input setup

% input folder: extracted raw rest data with grayMatter

inputfolder = fullfile(codecorresParentfolder, SKTDataSubfolder);

image_type = 'tif';

f_AOI = [8 40];


shuffleN_psedoTest = 100;

%% Code start here
cond_cell = cond_cell_extract(animal);

EventPhases = SKT_eventPhases_extract();
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal);

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
notAOI_chns = notInterested_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns

for ei = 1: length(EventPhases)
    event = EventPhases{ei};
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        
        subpdsavefolder = fullfile(savefolder, pdcond);
        if ~exist(subpdsavefolder, 'dir')
            mkdir(subpdsavefolder);
        end
        
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        % load(and extract) ciCohPhasefile
        ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond '_' event '_align2' align2name '.mat']);
        if(~exist(ciCohPhasefile, 'file'))
            
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end
            
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            if segVFolder
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            
            %  extract and save deltaphis_allChnsTrials and cicoh
            [deltaphis_allChnsTrials, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfptrials, fs, f_AOI);
            

            ntrials = size(lfptrials, 3);
            save(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
                      
            %  extract and save psedociCohs
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile);
            
            
            clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
            clear t_minmax_reach files lfptrials
        end
        load(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
        if ~exist('psedociCohs','var')
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile);
            load(ciCohPhasefile, 'psedociCohs');
            clear t_minmax_reach lfptrials
        end
        
        nshuffle = size(psedociCohs, 4);
        if(nshuffle < shuffleN_psedoTest)
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile);
            load(ciCohPhasefile, 'psedociCohs');
            clear t_minmax_reach lfptrials
        end
        nshuffle = size(psedociCohs, 4);
        
        % extract sigciCoh
        [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh);
       

        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [ciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(ciCoh, deltaphis_allChnsTrials, T_chnsarea);
        [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, ciCoh_flatten, deltaphis_flatten, removed_chns);
        
        
        % plot and save ciCoh Histogram image
        titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
        plot_ciCohHistogram(ciCoh_flatten_used, chnPairNames_used, f_selected, titlename);
        saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        clear titlename  saveimgname
        
        
        
        % rose histogram of deltaphis_allChnsTrials
        titlename_prefix = [animal '-'  pdcond '-'  event];
        subtitlename = [event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s, align2 = ' char(align2) ', ntrials = ' num2str(ntrials)];
        savefile_prefix = [animal 'trialPhaseDiff'];
        savefile_suffix = [event '_' pdcond '_align2' char(align2)];
        plotsave_deltaphirose(deltaphis_flatten_used, ciCoh_flatten_used, chnPairNames_used, f_selected, titlename_prefix, subtitlename, subpdsavefolder, savefile_prefix, savefile_suffix, image_type);
        clear titlename_prefix subtitlename savefile_prefix savefile_suffix
        
         
        clear pdcond subpdsavefolder align2 t_AOI align2name ciCohPhasefile
        clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
    end
end

function [deltaphis_allChnsTrials, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfptrials, fs, f_AOI)
%
%
%   
%   Output:
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials
%       ciCoh:  nchns * nchns * nf
%       f_selected: nf * 1



% extract ciCoh
[ciCoh, f_selected] = ciCohSKTAllchns_FFT_NoAmp(lfptrials, fs, f_AOI);

% extract deltaphis for all chns and trials
[nchns, ~, ntrials] = size(lfptrials);
nf = length(f_selected);
deltaphis_allChnsTrials = zeros(nchns, nchns, nf, ntrials);
for chni = 1: nchns - 1
    lfptriali = squeeze(lfptrials(chni, :, :));
    [phisi, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfptriali, fs, f_AOI);
    for chnj = chni + 1 : nchns
        lfptrialj = squeeze(lfptrials(chnj, :, :));
        [phisj, ~, ~] = phaseAmp_SKTPerTrial_FFT(lfptrialj, fs, f_AOI);
        deltaphis_allChnsTrials(chni, chnj, :, :) = phisi - phisj;
        clear lfptrialj phisj
    end
    clear lfptriali phisi
end


function [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh)
%
%   1. fit a normal distribution based on psedociCohs
%
%   2. extract pvalues using permutation test and extract significate ciCoh using Benjamini & Hochberg (1995) procedure 
%
%   Input:
%       psedociCohs: nchns * nchns * nf * nshuffles
%
%   Output:
%       sigciCoh: nchns * nchns * nf (if not sig, set 0; otherwise remain)

% fit a normal distribution
[nchns, ~, nf, ~] = size(psedociCohs);
mus = zeros(nchns, nchns, nf);
stds = zeros(nchns, nchns, nf);
for fi = 1: nf
    for chni = 1 : nchns -1
        for chnj = chni + 1 : nchns
            x = squeeze(psedociCohs(chni, chnj, nf, :));
            pd = fitdist(x,'Normal');
            mus(chni, chnj, fi) = pd.mu;
            stds(chni, chnj, fi) = pd.std;
            clear x pd
        end
    end
end

% pvalues using permutation test
[nchns, ~, nf] = size(ciCoh);
pvals = zeros(size(ciCoh));
for fi = 1 : nf
    for chni = 1: nchns -1
        for chnj = chni : nchns
            mu1 = mus(chni, chnj, fi);
            std1 = stds(chni, chnj, fi);
            pd = makedist('Normal','mu',mu1,'sigma',std1);
            
            x = ciCoh(chni, chnj, fi);
            pvals(chni, chnj, fi) = (1 - cdf(pd,abs(x)));
            
            clear x
            clear mu1 std1 pd
        end
    end
end
% Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR)
[h, ~, ~, ~]=fdr_bh(pvals);

% set values not significant as 0
sigciCoh = ciCoh;
sigciCoh(h==0) = 0;



function [ciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(ciCoh, deltaphis_allChnsTrials, T_chnsarea)
%
%   Input:
%       ciCoh: nchns * nchns * nf
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials
%       T_chnsarea: nchns * 5
%
%
[nchns, ~, nf] = size(ciCoh);
chnPairNames = {};
ciCoh_flatten = zeros(nchns * (nchns -1)/2, nf);
ntrials = size(deltaphis_allChnsTrials, 4);
deltaphis_flatten = zeros(nchns * (nchns -1)/2, nf, ntrials);
ci = 0;
for chni = 1 : nchns -1
    for chnj = chni + 1  : nchns
        chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
        
        ci = ci + 1;
        ciCoh_flatten(ci, :) = ciCoh(chni, chnj, :);
        deltaphis_flatten(ci, :, :) = deltaphis_allChnsTrials(chni, chnj, :, :);
    end
end


function [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, ciCoh_flatten, deltaphis_flatten, removed_chns)
%
% Input:
%   ciCoh_flatten: npairs * nf
%   deltaphis_flatten: npairs * nf * ntrials
%   

removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
clear M1DBS_mask STN2GP_mask

ciCoh_flatten_used = ciCoh_flatten(usedChnPairsMask, :);
deltaphis_flatten_used = deltaphis_flatten(usedChnPairsMask, :, :);
chnPairNames_used = chnPairNames(usedChnPairsMask);



function plot_ciCohHistogram(ciCoh_flatten, chnPairNames, f_selected, titlename)
fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;


% plot
figure;
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
imagesc(ciCoh_flatten)
colormap(jet)
set(gca, 'Position', [0.09 0.05 0.9 0.88])
[npairs, nf] = size(ciCoh_flatten);
xticks([1:nf])
xticklabels(round(f_selected,2))
yticks([1:npairs]);
set(gca,'YTickLabel',chnPairNames,'fontsize',12,'FontWeight','bold')
xlabel('freqs')
title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
set(gca,'CLim', [0 1])
colorbar



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
    
    clear s_stn s_gp chnPair
end



function plotsave_deltaphirose(deltaphis_flatten, ciCoh_flatten, chnPairNames, f_selected, titlename_prefix, subtitlename, savefolder, savefile_prefix, savefile_suffix, image_type)
nbins = 10;

fig_left = 50;
fig_bottom = 50;
fig_width = 800;
fig_height = 800;

nchnPairs = length(chnPairNames);
nf = length(f_selected);
for chnPairi = 1 : nchnPairs
    chnPairName = chnPairNames{chnPairi};
    for nfi = 1 : nf
        deltaphi = squeeze(deltaphis_flatten(chnPairi, nfi, :));
        f = f_selected(nfi);
        
        figure;
        set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
        set(gca, 'Position', [0.05 0.05 0.85 0.85])
        polarhistogram(deltaphi, nbins);
        
        titlename = [titlename_prefix ' Trial Phase Diff of ' chnPairName ' at ' num2str(round(f)) 'Hz'];
        title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
        
        % subtitle
        annotation(gcf,'textbox',[0.5 0.017 0.5 0.032], 'String',{subtitlename}, 'LineStyle','none', 'FitBoxToText','off');
        
        % plot icoh if sig
        sig = false;
        icoh = ciCoh_flatten(chnPairi, nfi);
        if(icoh > 0)
            sig = true;
        end
        if sig
            text(0.8, 0.9, ['icoh = ' num2str(round(icoh, 3)) '*'], 'FontSize', 12);
        end
        
        % save
        sigstr = '';
        if sig
            sigstr = '_sig';
        end
        savefile =  fullfile(savefolder, [savefile_prefix '_pair' chnPairName '_' num2str(round(f))  'Hz_' savefile_suffix sigstr '.' image_type]);
        saveas(gcf,savefile, image_type);
        
        clear icoh
        close all
    end
end

function psedociCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohPhasefile)
%
% lfptrials: nchns * ntemp * ntrials
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile

nchns = size(lfptrials, 1);

load(ciCohPhasefile, 'psedociCohs');
if(~exist('psedociCohs', 'var'))
    psedociCohs = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs, 4) + 1;
end
disp(['psedo test start at ' num2str(shuffi_str) ' times'])
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
    psedociCohs = cat(4, psedociCohs, psedociCohs1Time);
    clear psedociCohs1Time
    if(mod(si,10) == 0)
        disp(['psedo test finished ' num2str(si) ' times'])
        save(ciCohPhasefile, 'psedociCohs', '-append');
    end
end

