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
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '/', '[A-Za-z]*']);
elseif ispc
    % Code to run on Windows platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '\\', '[A-Za-z]*']);
else
    disp('Platform not supported')
end
animal = codecorresfolder(fi + length('NHPs') + 1:j);

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

EventPhases = {'preMove'; 'Anticipated';'earlyReach'; 'Return';'lateReach'};

twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];


if strcmpi(animal, 'bug')
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
end
if strcmpi(animal, 'jo')
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
    
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

[t_minmax_reach_normal, ~, t_minmax_reach_mild, ~, t_minmax_reach_moderate, ~] = goodSKTTrials_reachReturn_tcritiria(animal);
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    % each event Phase
    for ei = 1: length(EventPhases)
        event = EventPhases{ei};
        [align2, t_AOI, align2name] = SKTEventPhase_align2_tAOI_extract(event, animal, pdcond);
        
        
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
        
        [nchns, ntemp, ntrials] = size(lfptrials);
        
        % Rest data
        if ei == 1
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
        end
        
        [iCoh_trial, f_selected_trial] = imCohSKT_FFT_NormalizedAMP(lfptrials, twin, toverlap, fs_SKT, f_AOI, t_AOI, tdur_trial);
        
        % iCohChanges
        if any(f_selected_trial ~= f_selected_rest)
            disp(['f_selected_trial ~= f_selected_rest ' pdcond ' ' event])
            continue;
        end
        if ~isequal(T_chnsarea_rest, T_chnsarea_SKT)
            disp(['T_chnsarea_rest ~= T_chnsarea_SKT ' pdcond ' ' event ])
            continue;
        end
        iCohChanges_trial = iCoh_trial - iCoh_rest;
        
        % psedo Test
        shuffleN = 500;
        lfp_combined = cat(lfpdata_rest, lfptrials, 3);
        ntotal = size(lfp_combined, 3);
        for si = 1 : shuffleN
            randomSKTInds =  randsample(ntotal,ntrials);
            randomRestInds= arrayfun(@(x) any(randomSKTInds==x),[1: ntotal]);
            psedolfp_SKT = lfp_combined(:, :, randomSKTInds);
            psedolfp_Rest = lfp_combined(:, :, randomRestInds);
            [psedoiCoh_SKT, ~] = imCohSKT_FFT_NormalizedAMP(psedolfp_SKT, twin, toverlap, fs_SKT, f_AOI, t_AOI, tdur_trial);
            [psedoiCoh_rest, ~] = imCohRest_FFT_NormalizedAMP(psedolfp_Rest, twin, toverlap, fs_rest, f_AOI);
            
            if ~exist('psedoiCohChange', 'var')
                psedoiCohChanges = psedoiCoh_SKT - psedoiCoh_rest;
            else
                psedoiCohChanges = cat(4, psedoiCohChanges, psedoiCoh_SKT - psedoiCoh_rest);
            end
            
            clear randomSKTInds randomRestInds psedolfp_SKT psedoiCoh_rest
            clear psedoiCoh_SKT psedoiCoh_rest
        end
        
        % fit a normal distribution to psedoiCohChanges for each chni-chnj pair
        nf = size(psedoiCohChanges, 3);
        mus = zeros(chns, chns, nf);
        stds = zeros(chns, chns, nf);
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
        [nchns, ~, nf] = size(iCohChanges_trial);
        pvals = zeros(size(iCoh_acrossTime));
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
        [h, ~, ~, ~]=fdr_bh(pvals);
        
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
