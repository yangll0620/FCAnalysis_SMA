function m5_radarGraphofPhase()
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

presentName = 'trialPhaseDiff';

%%  input setup
inputfolder_SKT = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_goodReach');
% input folder: extracted raw rest data with grayMatter
% switch lower(animal)
%     case 'jo'
%         inputfolder_SKT = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');
%     case 'kitty'
%         inputfolder_SKT = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_goodReach');
% end

[pathstr,~,~] = fileparts( codecorresParentfolder );
inputfolder_Rest = fullfile(pathstr, 'Rest', 'm3_restData_rmChns_avgArea');

EventPhases = {'preMove'; 'Anticipated';'earlyReach'; 'lateReach'};

twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];


fig_left = 50;
fig_bottom = 50;
fig_width = 600;
fig_height = 600;

image_type = 'tif';

% plot setup
ylimvec = [-pi pi];
ytickvec = [-pi: pi/2 : pi];
yticklabelvec = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
xlabelname = 'cos \Delta\Phi';
ylabelname = 'sin \Delta\Phi';

color_separate = 'k';

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
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
end


%% Code Start Here
cond_cell = cond_cell_extract(animal);

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns];
clear unwanted_DBS noisy_chns


tic
[t_minmax_reach_normal, ~, t_minmax_reach_mild, ~, t_minmax_reach_moderate, ~] = goodSKTTrials_reachReturn_tcritiria(animal);
for ci = 2 : length(cond_cell)
    pdcond = cond_cell{ci};
    subsavefolder = fullfile(savefolder, pdcond);
    if ~exist(subsavefolder, 'dir')
        mkdir(subsavefolder);
    end
    
    % each event Phase
    for ei = 1: length(EventPhases)
        event = EventPhases{ei};
        disp([animal '-' pdcond '-' event])
        [align2, t_AOI, align2name] = SKTEventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        savefile = fullfile(savefolder, [animal  '_' pdcond '_' event '_align2' char(align2) '.mat']);
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
            [lfptrials, fs_SKT, T_chnsarea_SKT] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, t_minmax_reach);
            
            [nchns, ~, ntrials] = size(lfptrials);
           
            
            for chni = 1: nchns -1
                lfptrialsi = squeeze(lfptrials(chni, :, :)); 
                for chnj = chni + 1 : nchns
                    lfptrialsj = squeeze(lfptrials(chnj, :, :));
                    [deltaPhis, deltaAmps, f_selected] = deltaPhi_eachTrial(lfptrialsi, lfptrialsj, fs_SKT, twin, toverlap, f_AOI, t_AOI, tdur_trial);
                    
                    if ~exist('deltaPhis_allchns', 'var')
                        nfs = size(deltaPhis, 1);
                        deltaPhis_allchns = zeros(nchns, nchns, nfs, ntrials);
                        deltaAmps_allchns = zeros(nchns, nchns, nfs, ntrials);
                        clear nfs
                    end
                    deltaPhis_allchns(chni, chnj, :, :) = deltaPhis;
                    deltaAmps_allchns(chni, chnj, :, :) = deltaAmps;
                    clear deltaPhis deltaAmps lfptrialsj
                end
                clear lfptrialsi
            end
            
            
            % show and save iCoh images
            % generate chnPairNames, such as M1-stn0-1
            chnPairNames = {};
            nfs = size(deltaPhis_allchns, 3);
            deltaPhis_trialsFlat = zeros(nchns * (nchns -1)/2, nfs, ntrials);
            deltaAmps_trialsFlat = zeros(nchns * (nchns -1)/2, nfs, ntrials);
            ci = 0;
            for chni = 1 : nchns -1
                for chnj = chni + 1  : nchns
                    chnPairNames = [chnPairNames; {[T_chnsarea_SKT.brainarea{chni} '-'  T_chnsarea_SKT.brainarea{chnj}]}];
                    
                    ci = ci + 1;
                    deltaPhis_trialsFlat(ci, :, :) = deltaPhis_allchns(chni, chnj, :, :);
                    deltaAmps_trialsFlat(ci, :, :) = deltaAmps_allchns(chni, chnj, :, :);
                end
            end
            clear nfs ci deltaPhis_allchns deltaAmps_allchns
            
            
            % save data
            save(savefile, 'deltaPhis_trialsFlat', 'deltaAmps_trialsFlat', 'f_selected',  'chnPairNames')
            
            clear t_minmax_reach tdur_trial
            clear lfptrials fs_SKT T_chnsarea_SKT iCoh_trial f_selected_trial
            clear lfpdata_rest fs_rest T_chnsarea_rest iCoh_rest f_selected_rest 
            clear nchns ntrials nf ntotal
            clear lfp_combined psedoiCohChanges mus stds pvals h
            clear ci iCohChanges_trial
            clear iCohChanges_trialsFlat f_selected  T_chnsarea chnPairNames 
        end
        
        % deltaPhis_trialsFlat: nchnPairs * nfs * ntrials
        load(savefile, 'deltaPhis_trialsFlat', 'deltaAmps_trialsFlat', 'f_selected', 'chnPairNames');
        
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
        clear M1DBS_mask STN2GP_mask
        
        showData = deltaPhis_trialsFlat(usedChnPairsMask, :, :); 
        chnPairNames_show = chnPairNames(usedChnPairsMask);
        
        % plot & save
        [nchnPairs, nfs, ntrials] = size(showData); % showData: nchnPairs * nfs * ntrials
        for chnPairi = 1 : nchnPairs
            chnPairName = chnPairNames_show{chnPairi};
            for nfi = 1 : nfs
                phitmp = squeeze(showData(chnPairi, nfi, :));
                x = cos(phitmp);
                y = sin(phitmp);
                f = f_selected(nfi);
                
                figure;
                set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
                set(gca, 'Position', [0.05 0.05 0.85 0.85])
                               
                plot(x, y, '.')
                set(gca, 'XLim', [-1 1], 'YLim', [-1 1]);
                set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
                xlabel(xlabelname)
                ylabel(ylabelname)
                

                titlename = [animal '-'  pdcond '-'  event ' Trial Phase Differences of ' ...
                    chnPairName ' at ' num2str(round(f)) 'Hz'];
                title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
                
                % subtitle
                titleH = get(gca, 'title');
                titlePos = get(titleH, 'Position');
                set(titleH, 'Position', [round(titlePos(1)), titlePos(2)+ 0.3]); 
                subtitlename = [event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s, align2 = ' char(align2) ', ntrials = ' num2str(ntrials)];
                subtitlepos = [round(titlePos(1)), titlePos(2)+ 0.1];
                text(subtitlepos(1), subtitlepos(2), subtitlename, 'FontSize', 12);
                
                % save
                savefile =  fullfile(subsavefolder, [animal presentName '_pair' chnPairName '_' num2str(round(f))  'Hz_' event '_' pdcond '_align2' char(align2) '.' image_type]);
                saveas(gcf,savefile, image_type);
                
                close all
                
                clear x f titlename 
                clear titleH titlePos subtitlename subtitlepos
                clear savefile
               
            end
            
            clear chnPairName
        end

        
    end
    
    clear subsavefolder
end


function [deltaPhis, deltaAmps, f_selected] = deltaPhi_eachTrial(lfptrialsi, lfptrialsj, fs, twin, toverlap, f_AOI, t_AOI, tdur_trial)
%
%   Input:
%       lfptrialsi, lfptrialsj: lfp data for channel i and j (ntemp * ntrials)
%
%       f_AOI, t_AOI: frequencies and time duration of interest, i.e. f_AOI
%       = [8 40], t_AOI = [-0.2 0]
%
%       twin, toverlap: time duration to calculate ciCOH, i.e. twin = 0.2; toverlap = 0.15;
%
%       tdur_trial: time duration for used trial data, i.e tdur_trial = [-0.5 0.5]
%
%   Outputs:
%       deltaPhis: deltaPhi for each trial (nfs  * ntrials)

deltaPhis = [];
deltaAmps = [];
[~, ntrials] = size(lfptrialsi);
nwin = twin * fs;
noverlap = (toverlap * fs);
for triali = 1: ntrials
    x = lfptrialsi(: , triali);
    y = lfptrialsj(: , triali);
    
    [Sx, fx, tx, ~] = spectrogram(x, nwin, noverlap,[],fs); % Sx: nf * nt
    [Sy, ~, ~, ~] = spectrogram(y, nwin, noverlap,[],fs); % Sy: nf * nt
    
    if triali == 1
        freqs = fx;
        times = tx + tdur_trial(1);
        idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
        idx_t = (times>=t_AOI(1) &  times<=t_AOI(2));
        f_selected = round(freqs(idx_f), 2);
        clear freqs times
    end
    
    phix = angle(Sx(idx_f, idx_t));
    phiy = angle(Sy(idx_f, idx_t));
    ampx = abs(Sx(idx_f, idx_t));
    ampy = abs(Sy(idx_f, idx_t));
    deltaPhi = mean(phix - phiy, 2);
    deltaAmp = mean(ampx ./ ampy, 2);
    deltaPhis = cat(2, deltaPhis, deltaPhi);
    deltaAmps = cat(2, deltaAmps, deltaAmp);
    
    
    clear x y Sx fx tx Sy
    clear phix phiy deltaPhi
    clear ampx ampy deltaAmp
end