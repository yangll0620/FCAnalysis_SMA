function m5_imCohUsingFFT_EventPhase()
%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);

%% save setup
savefolder = codecorresfolder;


%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');


fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;

EventPhases = {'preMove'; 'Anticipated';'earlyReach'; 'Return';'lateReach'};

image_type = 'tif';


twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];


%% Code start here
cond_cell = cond_cell_extract(animal);
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);
if strcmpi(animal, 'bug')
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
end
if strcmpi(animal, 'jo')
    tdur_trial_normal = [-1 0.5];
    tdur_trial_mild = [-1 0.5];
    tdur_trial_moderate = [-1 0.5];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    tdur_trial_normal = [-1 0.5];
    tdur_trial_moderate = [-1 0.5];
end

if strcmpi(animal, 'pinky')
    tdur_trial_normal = [-1 0.5];
    tdur_trial_mild = [-1 0.5];
    tdur_trial_moderate = [-1 0.5];
end



unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
notInterested_chns = notInterested_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns notInterested_chns];
clear unwanted_DBS noisy_chns notInterested_chns

for ei = 1: length(EventPhases)
    event = EventPhases{ei};
    
    
    % SKT phase
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        [align2, t_AOI, align2name] = SKTEventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        
        savefile_prefix = fullfile(savefolder, [animal '_' pdcond ]);
        savefile_SKT = [savefile_prefix '_' event '_align2' align2name '.mat'];

        if ~exist(savefile_SKT)
            if strcmp(pdcond, 'normal')
                t_minmax_reach = t_minmax_reach_normal;
                t_minmax_return = t_minmax_return_normal;
                tdur_trial = tdur_trial_normal;
            else
                if strcmp(pdcond, 'mild')
                    t_minmax_reach = t_minmax_reach_mild;
                    t_minmax_return = t_minmax_return_mild;
                    tdur_trial = tdur_trial_mild;
                else
                    if strcmp(pdcond, 'moderate')
                        t_minmax_reach = t_minmax_reach_moderate;
                        t_minmax_return = t_minmax_return_moderate;
                        tdur_trial = tdur_trial_moderate;
                    end
                end
                
            end
            
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end
            
            [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach);
            
            [~, ~, ntrials] = size(lfptrials);
            
            % calculate imaginary of Coherency
            [iCoh_acrossTime, f_selected] = imCohSKT_FFT_NormalizedAMP(lfptrials, twin, toverlap, fs, f_AOI, t_AOI, tdur_trial);
            
            lfp1 = squeeze(lfptrials(1, :, :));
            lfp2 = squeeze(lfptrials(10, :, :));
            [mus, stds] = psedo_SKTLFP_Test_acrossTime(lfp1, lfp2, 100, fs, twin, toverlap, f_AOI, t_AOI - tdur_trial(1));
            clear lfp1 lfp2
            
            
            % pvalues using permutation test
            [nchns, ~, nf] = size(iCoh_acrossTime);
            pvals = zeros(size(iCoh_acrossTime));
            for fi = 1 : nf
                
                mu = mus(fi,1);
                std = stds(fi, 1);
                pd = makedist('Normal','mu',mu,'sigma',std);
                
                for chni = 1: nchns -1
                    for chnj = chni : nchns
                        x = iCoh_acrossTime(chni, chnj, fi);
                        pvals(chni, chnj, fi) = (1 - cdf(pd,abs(x)));
                        pvals(chnj, chni, fi) = pvals(chni, chnj, fi);
                        clear x
                    end
                end
                
                clear mu std pd
            end
            % Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR)
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
            
            % set values not significant as 0
            iCoh_acrossTime(h==0) = 0;
            
            % show and save iCoh images
            
            % generate chnPairNames, such as M1-stn0-1
            nf = size(iCoh_acrossTime, 3);
            chnPairNames = {};
            iCoh_1time = zeros(nchns * (nchns -1)/2, nf);
            ci = 0;
            for chni = 1 : nchns -1
                for chnj = chni + 1  : nchns
                    chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
                    
                    ci = ci + 1;
                    iCoh_1time(ci, :) = iCoh_acrossTime(chni, chnj, :);
                end
            end
            
            % save data
            save(savefile_SKT, 'iCoh_1time', 'f_selected', 'chnPairNames', 'ntrials')
          
        else
            load(savefile_SKT, 'iCoh_1time', 'f_selected',  'chnPairNames', 'ntrials')
        end
        
        
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames); 
        usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
        clear M1DBS_mask STN2GP_mask
        
        showData = abs(iCoh_1time(usedChnPairsMask, :));
        chnPairNames_show = chnPairNames(usedChnPairsMask);
        
        % plot
        figure;
        set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
        imagesc(showData)
        colormap(jet)
        set(gca, 'Position', [0.09 0.05 0.9 0.88])
        [npairs, nf] = size(showData);
        xticks([1:nf])
        xticklabels(round(f_selected,2))
        yticks([1:npairs]);
        set(gca,'YTickLabel',chnPairNames_show,'fontsize',12,'FontWeight','bold')
        xlabel('freqs')
        title([animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' 'align2 = ' align2name ', ntrials = ' num2str(ntrials)], ...
            'FontSize', 15, 'FontWeight', 'normal')
        set(gca,'CLim', [0 1])
        colorbar
        
        
        
        chnPair_prev = '';
        for ci = 1: length(chnPairNames_show)
            chnPair = chnPairNames_show{ci};
            
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
        
        % save image
        saveas(gcf, fullfile(savefolder, [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type]), image_type);
       
        close all
        
        clear pdcond t_minmax_reach t_minmax_return tdur_trial
        clear files lfptrials fs T_chnsarea nchns ntrials
        clear iCoh iCoh_acrossTime mus stds
        clear idx_f idx_t f_selected t_selected
        clear nf chnPairNames iCoh_1time ci
        clear M1DBS_mask STN2GP_mask usedChnPairsMask showData npairs
        clear savefile_prefix
    end
    
    
    clear t_AOI align2 event
end

%% copy code to savefolder
[~, codefilename]= fileparts(codefilepath);
status = copyfile([codefilepath '.m'], fullfile(savefolder, [codefilename '.m']));
disp(['copied code status = ' num2str(status)])
