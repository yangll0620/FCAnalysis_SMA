function m4_imCohChangesUsingFFT()
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
inputfolder_SKT = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_goodReach');
files = dir(fullfile(inputfolder_SKT, ['*.mat']));
if(isempty(files))
    inputfolder_SKT = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');
end


[pathstr,~,~] = fileparts( codecorresParentfolder );
inputfolder_Rest = fullfile(pathstr, 'Rest', 'm3_restData_rmChns_avgArea');

EventPhases = {'preMove'; 'Anticipated';'earlyReach'; 'lateReach'};

twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];


fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;

image_type = 'tif';

shuffleN_psedoTest = 500;

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


%% Code Start Here
cond_cell = cond_cell_extract(animal);

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns];
clear unwanted_DBS noisy_chns


tic
[t_minmax_reach_normal, ~, t_minmax_reach_mild, ~, t_minmax_reach_moderate, ~] = goodSKTTrials_reachReturn_tcritiria(animal);
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    % each event Phase
    for ei = 1: length(EventPhases)
        event = EventPhases{ei};
        disp([animal '-' pdcond '-' event])
        [align2, t_AOI, align2name] = SKTEventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        savefile = fullfile(savefolder, [animal '_shuffle' num2str(shuffleN_psedoTest) '_' pdcond '_' event '_align2' char(align2) '.mat']);
        if ~exist(savefile)
            switch lower(pdcond)
                case 'normal'
                    t_minmax_reach = t_minmax_reach_normal;
                    tdur_trial = tdur_trial_normal;
                case 'mild'
                    t_minmax_reach = t_minmax_reach_mild;
                    tdur_trial = tdur_trial_mild;
                case 'moderate'
                    t_minmax_reach = t_minmax_reach_moderate;
                    tdur_trial = tdur_trial_moderate;
            end
            
            files = dir(fullfile(inputfolder_SKT, ['*_' pdcond '_*.mat']));       
            switch lower(animal)
                case 'pinky'
                    [lfptrials, fs_SKT, T_chnsarea_SKT] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach);
                otherwise
                    [lfptrials, fs_SKT, T_chnsarea_SKT] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, t_minmax_reach);
            end
                   
            
            [nchns, ~, ntrials] = size(lfptrials);
            
            if ~exist('iCoh_rest', 'var') % deal with  Rest data once
                disp(['iCoh_rest Calculate in ' pdcond]);
                files_Rest = dir(fullfile(inputfolder_Rest, ['*_' pdcond '_*.mat']));
                [lfpdata_rest, fs_rest, T_chnsarea_rest]= seg2ShortSegments(files_Rest, tdur_trial(2) - tdur_trial(1));
                
                if fs_SKT ~= fs_rest
                    disp(['fs_SKT ~= fs_rest ' pdcond])
                    continue;
                end
                
                % match the trial number in rest and SKT
                nsegs = size(lfpdata_rest, 3);
                randomSKTInds =  randsample(nsegs,ntrials);
                lfpdata_rest = lfpdata_rest(:, :, randomSKTInds);
                
                [iCoh_rest, f_selected_rest] = imCohRest_FFT_NormalizedAMP(lfpdata_rest, twin, toverlap, fs_rest, f_AOI);
                clear files_Rest nsegs randomInds
                toc
            end
            
            [iCoh_trial, f_selected_trial] = imCohSKT_FFT_NormalizedAMP(lfptrials, twin, toverlap, fs_SKT, f_AOI, t_AOI, tdur_trial);
            
            % check consistent between trial and rest for f_selected and T_chnsarea
            if any(f_selected_trial ~= f_selected_rest)
                disp(['f_selected_trial ~= f_selected_rest ' pdcond ' ' event])
                continue;
            end
            if ~isequal(T_chnsarea_rest.brainarea, T_chnsarea_SKT.brainarea)
                disp(['T_chnsarea_rest ~= T_chnsarea_SKT ' pdcond ' ' event ])
                continue;
            end
            f_selected = f_selected_trial;
            T_chnsarea = T_chnsarea_SKT;
            
            % iCohChanges_trials
            iCohChanges_trial = iCoh_trial - iCoh_rest;
            
            nf = size(iCohChanges_trial, 3);
            
            % psedo Test
            lfp_combined = cat(3, lfpdata_rest, lfptrials);
            ntotal = size(lfp_combined, 3);
            for si = 1 : shuffleN_psedoTest
                if(mod(si, 10) == 0)
                    disp(['pesdo test ' num2str(si)])
                    toc
                end
                randomSKTInds =  randsample(ntotal,ntrials);
                randomRestInds= arrayfun(@(x) any(randomSKTInds==x),[1: ntotal]);
                psedolfp_SKT = lfp_combined(:, :, randomSKTInds);
                psedolfp_rest = lfp_combined(:, :, randomRestInds);
                [psedoiCoh_SKT, ~] = imCohSKT_FFT_NormalizedAMP(psedolfp_SKT, twin, toverlap, fs_SKT, f_AOI, t_AOI, tdur_trial);
                [psedoiCoh_rest, ~] = imCohRest_FFT_NormalizedAMP(psedolfp_rest, twin, toverlap, fs_rest, f_AOI);
                
                if ~exist('psedoiCohChanges', 'var')
                    psedoiCohChanges = psedoiCoh_SKT - psedoiCoh_rest;
                else
                    psedoiCohChanges = cat(4, psedoiCohChanges, psedoiCoh_SKT - psedoiCoh_rest);
                end
                
                clear randomSKTInds randomRestInds psedolfp_SKT psedoiCoh_rest
                clear psedoiCoh_SKT psedoiCoh_rest
            end
            
            % fit a normal distribution to psedoiCohChanges for each chni-chnj pair
            disp('fit a normal distribution')
            mus = zeros(nchns, nchns, nf);
            stds = zeros(nchns, nchns, nf);
            for chni = 1 : nchns-1
                for chnj = chni + 1 : nchns
                    for fi = 1 : nf
                        pd = fitdist(squeeze(psedoiCohChanges(chni, chnj, fi, :)),'Normal');
                        mus(chni, chnj, fi) = pd.mu;
                        stds(chni, chnj, fi) = pd.sigma;
                        
                        clear pd
                    end
                end
            end
            
            % pvalues using permutation test
            disp('pvalues using permutation test')
            pvals = zeros(size(iCohChanges_trial));
            for fi = 1 : nf
                for chni = 1: nchns -1
                    mu = mus(chni, chnj, fi);
                    std = stds(chni, chnj, fi);
                    pd = makedist('Normal','mu',mu,'sigma',std);
                    for chnj = chni : nchns
                        x = iCohChanges_trial(chni, chnj, fi);
                        pvals(chni, chnj, fi) = (1 - cdf(pd,abs(x)));
                        pvals(chnj, chni, fi) = pvals(chni, chnj, fi);
                        clear x
                    end
                    clear mu std pd
                end
            end
            % Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR)
            [h, ~, ~, ~] = fdr_bh(pvals);
            
            % set values not significant as 0
            iCohChanges_trial(h==0) = 0;
            
            
            % show and save iCoh images
            % generate chnPairNames, such as M1-stn0-1
            chnPairNames = {};
            iCohChanges_trialsFlat = zeros(nchns * (nchns -1)/2, nf);
            ci = 0;
            for chni = 1 : nchns -1
                for chnj = chni + 1  : nchns
                    chnPairNames = [chnPairNames; {[T_chnsarea_SKT.brainarea{chni} '-'  T_chnsarea_SKT.brainarea{chnj}]}];
                    
                    ci = ci + 1;
                    iCohChanges_trialsFlat(ci, :) = iCohChanges_trial(chni, chnj, :);
                end
            end
            
            % save data
            save(savefile, 'iCohChanges_trialsFlat', 'f_selected', 'chnPairNames', 'ntrials')
            
            clear t_minmax_reach tdur_trial
            clear lfptrials fs_SKT T_chnsarea_SKT iCoh_trial f_selected_trial
            clear lfpdata_rest fs_rest T_chnsarea_rest iCoh_rest f_selected_rest 
            clear nchns ntrials nf ntotal
            clear lfp_combined psedoiCohChanges mus stds pvals h
            clear ci iCohChanges_trial
            clear iCohChanges_trialsFlat f_selected  T_chnsarea chnPairNames 
        end
        
        load(savefile, 'iCohChanges_trialsFlat', 'f_selected',  'chnPairNames', 'ntrials')
        
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
        clear M1DBS_mask STN2GP_mask
        
        showData = iCohChanges_trialsFlat(usedChnPairsMask, :);
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
        title([animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s FC Changes,' 'align2 = ' char(align2) ', ntrials = ' num2str(ntrials)], ...
            'FontSize', 15, 'FontWeight', 'normal')
        set(gca,'CLim', [-0.5 0.5])
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
        saveas(gcf, fullfile(savefolder, [animal ' FCChanges_shuffle' num2str(shuffleN_psedoTest) '_' event '_' pdcond '_align2' char(align2) '.' image_type]), image_type);
        close all
    end
end


function [lfpdata, fs, T_chnsarea]= seg2ShortSegments(files, twin)
lfpdata = [];
for fi = 1: length(files)
    loadfilename = files(fi).name;
    load(fullfile(files(fi).folder, loadfilename), 'data_segments', 'fs', 'T_chnsarea');
    
    nwin = round(twin * fs);
    for segi = 1 : length(data_segments)
        seglfp = data_segments(segi).lfp;
        
        len = size(seglfp, 2);
        shortSegn = floor(len / nwin);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nwin + 1;
            endi = shortsegi * nwin;
            lfpdata = cat(3, lfpdata, seglfp(:, stri : endi));
            clear stri endi
        end
        clear seglfp len shortSegn shortsegi
    end
    
    if ~exist('fs_unit', 'var')
        fs_unit = fs;
    else
        if(fs_unit ~=fs)
            dis(['fs_unit ~=fs for file ' loadfilename])
        end
    end
    
    if ~exist('T_chnsarea_unit', 'var')
        T_chnsarea_unit = T_chnsarea;
    else
        if(~isequal(T_chnsarea_unit, T_chnsarea))
            dis(['T_chnsarea_unit ~= T_chnsarea for file ' loadfilename])
        end
    end
    
    clear nwin segi
    clear('data_segments', 'fs', 'T_chnsarea')
    clear loadfilename
    
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;