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
addpath(genpath(fullfile(codefolder,'toolbox')));

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

ciCohPhasefile_prefix =[animal ' ciCohPhasefile'];

%%  input setup

% input folder: extracted raw rest data with grayMatter

inputfolder = fullfile(codecorresParentfolder, SKTDataSubfolder);

image_type = 'tif';

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

for ei = 1: length(EventPhases)
    event = EventPhases{ei};
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        
        subpdsavefolder = fullfile(savefolder, pdcond);
        if ~exist(subpdsavefolder, 'dir')
            mkdir(subpdsavefolder);
        end
        
        disp([codefilename ' ' animal '-' event '-' pdcond])
        
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
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            
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
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            
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
            
            % remove unused chns
            removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
            lfptrials = lfptrials(~removedChns_mask, :, :);
            T_chnsarea = T_chnsarea(~removedChns_mask, :);
            
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile);
            load(ciCohPhasefile, 'psedociCohs');
            clear t_minmax_reach lfptrials
        end
        nshuffle = size(psedociCohs, 4);
        
        % extract sigciCoh
        [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh);
        

        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [sigciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCoh, deltaphis_allChnsTrials, T_chnsarea);
        [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, sigciCoh_flatten, deltaphis_flatten, removed_chns);
        
        
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
