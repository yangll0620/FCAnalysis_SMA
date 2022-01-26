function m6_imCohChangesUsingFFT_basedNormal_unifiedNHP(animal, varargin)


codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'toolbox')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));

% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
segVFolder = false;
if strcmpi(animal, 'Kitty')
    segVFolder = true;
end

if segVFolder
    SKTSubfolder = 'SKT_SegV';
    SKTlfpSubfolder = 'm2_segSKTData_SelectTrials_goodReach';
else
    SKTSubfolder = 'SKT';
    SKTlfpSubfolder = 'm2_SKTData_SelectTrials';
end
SKTCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, codefilename);
[SKTCodecorresfolder, SKTCodecorresParentfolder] = code_corresfolder(SKTCodefilepath, true, false);
pipelinefolder = fullfile(fileparts(codefolder), 'pipeline');

%% save setup
savefolder = SKTCodecorresfolder;
copyfile2folder(codefilepath, savefolder);

ciCohChangesfile_prefix =[animal ' ciCohChangesfile'];

%%  input setup
inputfolder_ciCohSKT = fullfile(SKTCodecorresParentfolder, 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');
inputfolder_lfpSKT = fullfile(SKTCodecorresParentfolder, SKTlfpSubfolder);


image_type = 'tif';

twin = 0.2;
f_AOI = [8 40];


shuffleN_psedoTest = 500;

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



% each event Phase
for ei = 3: length(EventPhases)
    event = EventPhases{ei};
    
    % base lfp is normal
    basepd = 'normal';
    [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, basepd);
    load(fullfile(inputfolder_ciCohSKT, [animal ' ciCohPhasefile_' basepd '_' event '_align2' align2name '.mat']), 'ciCoh');
    ciCoh_base = ciCoh;
    
    files = dir(fullfile(inputfolder_lfpSKT, ['*_' basepd '_*.mat']));
    eval(['t_minmax_reach = t_minmax_reach_' basepd ';']);
    if segVFolder
        [lfptrials_base, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
    else
        [lfptrials_base, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
    end
    clear files
    
    % remove unused chns
    removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
    lfptrials_base = lfptrials_base(~removedChns_mask, :, :);
    T_chnsarea = T_chnsarea(~removedChns_mask, :);
    
    for ci = 2: length(cond_cell)
        pdcond = cond_cell{ci};
        
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond);
        disp([codefilename ': ' animal '-' pdcond '-' event])
        
        % load(and extract) ciCohPhasefile
        ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix  '_' pdcond '_' event '_align2' align2name '.mat']);
        if (~exist(ciCohChangesfile, 'file'))
            load(fullfile(inputfolder_ciCohSKT, [animal ' ciCohPhasefile_' pdcond '_' event '_align2' align2name '.mat']), 'ciCoh', 'T_chnsarea', 'f_selected', 'ntrials');
            ciCoh_pd = ciCoh;
            clear ciCoh
            
            ciCohChanges = ciCoh_pd - ciCoh_base;
            
            % save before psedo Test
            save(ciCohChangesfile, 'ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected');
            clear('ciCohChanges', 'T_chnsarea', 'f_selected', 'ntrials');
            
            % psedo Test
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder_lfpSKT, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, ~, ~] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, ~, ~] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            
            % remove unused chns
            lfptrials = lfptrials(~removedChns_mask, :, :);
            
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile);
            
            clear ciCohChanges files lfptrials
        end
        load(ciCohChangesfile, 'ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected', 'psedoiCohChanges');
        if ~exist('psedoiCohChanges','var') || size(psedoiCohChanges, 4) < shuffleN_psedoTest
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder_lfpSKT, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            clear T_chnsarea_lfprest
            
            %  extract and save psedociCohs
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile);
            
            clear lfptrials fs_lfpSKT T_chnsarea_lfpSKT files t_minmax_reach
            
            load(ciCohChangesfile, 'psedoiCohChanges');
        end
        nshuffle = size(psedoiCohChanges, 4);
        
        % extract sigciCoh
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        
        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        [ciCohChanges_flatten_used, chnPairNames_used]= ciCoh_Used(chnPairNames, sigciCohChanges_flatten, removed_chns);
        
        
        % plot and save ciCohChanges Histogram image
        titlename = [animal ' FC Changes '  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
        plot_ciCohHistogram(ciCohChanges_flatten_used, chnPairNames_used, f_selected, titlename, [-0.5 0.5]);
        saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        clear titlename  saveimgname
        close all
        
        clear align2 t_AOI align2name ciCohChangesfile
        clear('ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected', 'psedoiCohChanges');
    end
    
end

