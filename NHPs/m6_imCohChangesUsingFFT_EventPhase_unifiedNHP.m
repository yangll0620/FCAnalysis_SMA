function m6_imCohChangesUsingFFT_EventPhase_unifiedNHP(varargin)
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
inputfolder_ciCohRest = fullfile(SKTCodecorresParentfolder, 'm5_imCohUsingFFT_RestEqualDurSegnum_unifiedNHP');
inputfolder_ciCohSKT = fullfile(SKTCodecorresParentfolder, 'm4_imCohUsingFFT_EventPhase_unifiedNHP');
inputfolder_lfpSKT = fullfile(SKTCodecorresParentfolder, SKTlfpSubfolder);
inputfolder_lfpRest = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'Rest', 'm3_restData_rmChns_avgArea');

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


for ci = 3 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    % load Rest ciCoh
    load(fullfile(inputfolder_ciCohRest, [animal '_Rest_ciCohPhasefile_' pdcond '.mat']), 'ciCoh', 'T_chnsarea', 'f_selected', 'nsegs');
    ciCoh_Rest = ciCoh;
    T_chnsarea_ciCohRest = T_chnsarea;
    f_selected_ciCohRest = f_selected;
    clear ciCoh T_chnsarea f_selected
    
    % load lfpdata
    files = dir(fullfile(inputfolder_lfpRest, ['*_' pdcond '_*.mat']));
    [lfpsegs_rest, fs_lfprest, T_chnsarea_lfprest]= seg2ShortSegments(files, twin);
    % remove unused chns
    removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea_lfprest.brainarea);
    lfpsegs_rest = lfpsegs_rest(~removedChns_mask, :, :);
    T_chnsarea_lfprest = T_chnsarea_lfprest(~removedChns_mask, :);
    clear files removedChns_mask
    
    % check T_chnsarea equal
    if ~isequal(T_chnsarea_ciCohRest, T_chnsarea_lfprest)
        disp(['~isequal(T_chnsarea_ciCohRest, T_chnsarea_lfprest)'])
        continue;
    end
    T_chnsarea_rest = T_chnsarea_lfprest;
    clear T_chnsarea_lfprest T_chnsarea_ciCohRest
    
    
    
    % each event Phase
    for ei = 1: length(EventPhases)
        event = EventPhases{ei};
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        disp([codefilename ': ' animal '-' pdcond '-' event])
        
        % load(and extract) ciCohPhasefile
        ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix  '_' pdcond '_' event '_align2' align2name '.mat']);
        if (~exist(ciCohChangesfile, 'file'))

            % load SKT ciCoh
            load(fullfile(inputfolder_ciCohSKT, [animal ' ciCohPhasefile_' pdcond '_' event '_align2' align2name '.mat']), 'ciCoh', 'T_chnsarea', 'f_selected', 'ntrials');
            T_chnsarea_ciCohSKT = T_chnsarea;
            f_selected_ciCohSKT = f_selected;
            ciCoh_SKT = ciCoh;
            clear ciCoh T_chnsarea f_selected
            
            % iCohChanges
            if ~isequal(T_chnsarea_ciCohSKT.brainarea, T_chnsarea_rest.brainarea) || ~isequal(f_selected_ciCohSKT, f_selected_ciCohRest) || ntrials ~= nsegs
                disp('T_chnsarea_ciCohSKT or f_selected_ciCohSKT not equal');
                continue;
            end
            T_chnsarea = T_chnsarea_ciCohSKT;
            f_selected = f_selected_ciCohSKT;
            ciCohChanges = ciCoh_SKT - ciCoh_Rest;
            
            % save before psedo Test
            save(ciCohChangesfile, 'ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected');
            clear('ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected');
            
            % psedo Test
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder_lfpSKT, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, fs_lfpSKT, T_chnsarea_lfpSKT] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, fs_lfpSKT, T_chnsarea_lfpSKT] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea_lfpSKT.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea_lfpSKT = T_chnsarea_lfpSKT(~removedChns_mask, :);
            clear removedChns_mask
            
            if ~isequal(T_chnsarea_lfpSKT.brainarea, T_chnsarea_rest.brainarea) || ~isequal(fs_lfpSKT, fs_lfprest)
                disp('~isequal(T_chnsarea_lfpSKT.brainarea, T_chnsarea_rest.brainarea) || ~isequal(fs_lfpSKT, fs_lfprest)');
                continue;
            end
            fs = fs_lfpSKT;
            
            %  extract and save psedociCohs
            nsegs_lfprest = size(lfpsegs_rest, 3);
            ntrials = size(lfptrials, 3);
            randomSKTInds =  randsample(nsegs_lfprest,ntrials);
            lfpEqRest = lfpsegs_rest(:, :, randomSKTInds);
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials, lfpEqRest, fs, f_AOI, ciCohChangesfile);
            clear nsegs_lfprest ntrials randomSKTInds lfpEqRest
            
            clear lfptrials fs_lfpSKT T_chnsarea_lfpSKT files t_minmax_reach
        end
        load(ciCohChangesfile, 'ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected', 'psedoiCohChanges');
        
        if ~exist('psedoiCohChanges','var') || size(psedoiCohChanges, 4) < shuffleN_psedoTest
            eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
            files = dir(fullfile(inputfolder_lfpSKT, ['*_' pdcond '_*.mat']));
            if segVFolder
                [lfptrials, fs_lfpSKT, T_chnsarea_lfpSKT] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
            else
                [lfptrials, fs_lfpSKT, T_chnsarea_lfpSKT] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
            end
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea_lfpSKT.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea_lfpSKT = T_chnsarea_lfpSKT(~removedChns_mask, :);
            clear T_chnsarea_lfprest
            
            if ~isequal(T_chnsarea_lfpSKT.brainarea, T_chnsarea_rest.brainarea) || ~isequal(fs_lfpSKT, fs_lfprest)
                disp('~isequal(T_chnsarea_lfpSKT.brainarea, T_chnsarea_rest.brainarea) || ~isequal(fs_lfpSKT, fs_lfprest)');
                continue;
            end
            fs = fs_lfpSKT;
            
            %  extract and save psedociCohs
            nsegs_lfprest = size(lfpsegs_rest, 3);
            ntrials = size(lfptrials, 3);
            randomSKTInds =  randsample(nsegs_lfprest,ntrials);
            lfpEqRest = lfpsegs_rest(:, :, randomSKTInds);
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials, lfpEqRest, fs, f_AOI, ciCohChangesfile);
            clear nsegs_lfprest ntrials randomSKTInds lfpEqRest
            
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
        
        clear event align2 t_AOI align2name ciCohChangesfile
        clear('ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected', 'psedoiCohChanges');
    end
    
    
    clear pdcond f_selected_ciCohRest lfpsegs_rest fs_lfprest T_chnsarea_rest
end
