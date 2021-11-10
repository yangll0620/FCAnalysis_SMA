function m4_radarGraphofPhase()
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
fig_width = 1200;
fig_height = 600;

image_type = 'tif';


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
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
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
                    [deltaPhis, f_selected] = deltaPhi_eachTrial(lfptrialsi, lfptrialsj, fs_SKT, twin, toverlap, f_AOI, t_AOI, tdur_trial);
                    
                    if ~exist('deltaPhis_allchns', 'var')
                        nfs = size(deltaPhis, 1);
                        deltaPhis_allchns = zeros(nchns, nchns, nfs, ntrials);
                        clear nfs
                    end
                    deltaPhis_allchns(chni, chnj, :, :) = deltaPhis;
                    clear deltaPhis lfptrialsj
                end
                clear lfptrialsi
            end
            
            
            % show and save iCoh images
            % generate chnPairNames, such as M1-stn0-1
            chnPairNames = {};
            nfs = size(deltaPhis_allchns, 3);
            deltaPhis_trialsFlat = zeros(nchns * (nchns -1)/2, nfs, ntrials);
            ci = 0;
            for chni = 1 : nchns -1
                for chnj = chni + 1  : nchns
                    chnPairNames = [chnPairNames; {[T_chnsarea_SKT.brainarea{chni} '-'  T_chnsarea_SKT.brainarea{chnj}]}];
                    
                    ci = ci + 1;
                    deltaPhis_trialsFlat(ci, :, :) = deltaPhis_allchns(chni, chnj, :, :);
                end
            end
            clear nfs ci deltaPhis_allchns
            
            
            % save data
            save(savefile, 'deltaPhis_trialsFlat', 'f_selected', 'T_chnsarea_SKT', 'chnPairNames')
            
            clear t_minmax_reach tdur_trial
            clear lfptrials fs_SKT T_chnsarea_SKT iCoh_trial f_selected_trial
            clear lfpdata_rest fs_rest T_chnsarea_rest iCoh_rest f_selected_rest 
            clear nchns ntrials nf ntotal
            clear lfp_combined psedoiCohChanges mus stds pvals h
            clear ci iCohChanges_trial
            clear iCohChanges_trialsFlat f_selected  T_chnsarea chnPairNames 
        end
        
        load(savefile, 'deltaPhis_trialsFlat', 'f_selected',  'T_chnsarea_SKT', 'chnPairNames');
        
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
        clear M1DBS_mask STN2GP_mask
        
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