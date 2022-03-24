function m5_imCohUsingFFT_RestEqualDurSegnum_unifiedNHP(animal, varargin)
% 
%   Input
%       animal
%

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
copyfile2folder(codefilepath, savefolder);

ciCohPhasefile_prefix =[animal '_Rest_ciCohPhasefile'];

%%  input setup

% inputfolder_Rest & inputcodefolder_SKT
inputcodefolder_Rest = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'Rest', 'm3_restData_rmChns_avgArea');
[inputfolder_Rest, ~] = code_corresfolder(inputcodefolder_Rest, false, false);

inputcodefolder_SKT = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, 'm4_imCohUsingFFT_EventPhase_unifiedNHP');
[inputfolder_SKT, ~] = code_corresfolder(inputcodefolder_SKT, false, false);
clear inputcodefolder_Rest SKTSubfolder inputcodefolder_SKT

twin = 0.2;

image_type = 'tif';

f_AOI = [8 40];


shuffleN_psedoTest = 500;
runRose = false;

%% Code start here
cond_cell = cond_cell_extract(animal);

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
notAOI_chns = notInterested_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns

for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    subpdsavefolder = fullfile(savefolder, pdcond);
    if ~exist(subpdsavefolder, 'dir')
        mkdir(subpdsavefolder);
    end
    
    disp([codefilename ': ' animal '-' pdcond])
    
    % load(and extract) ciCohPhasefile
    ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond  '.mat']);
    if(~exist(ciCohPhasefile, 'file'))
        
        files = dir(fullfile(inputfolder_Rest, ['*_' pdcond '_*.mat']));
        [lfpsegs, fs, T_chnsarea]= seg2ShortSegments(files, twin);
        % remove unused chns
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(~removedChns_mask, :, :);
        T_chnsarea = T_chnsarea(~removedChns_mask, :);
        clear files removedChns_mask
        
        % match the trial number in rest and SKT
        load(fullfile(inputfolder_SKT, [animal ' ciCohPhasefile_' pdcond '_earlyReach_align2ReachOnset.mat']), 'ntrials')
        nsegs = size(lfpsegs, 3);
        randomSKTInds =  randsample(nsegs,ntrials);
        lfpsegs = lfpsegs(:, :, randomSKTInds);
        clear ntrials nsegs randomSKTInds
        
        
        
        %  extract and save deltaphis_allChnsTrials and cicoh
        [deltaphis_allChnsSegs, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
        
        % save
        nsegs = size(lfpsegs, 3);
        save(ciCohPhasefile, 'deltaphis_allChnsSegs', 'ciCoh', 'T_chnsarea', 'nsegs', 'f_selected');
            
        %  extract and save psedociCohs
        psedociCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile);
        
        clear files lfpsegs fs T_chnsarea
        clear('deltaphis_allChnsSegs', 'ciCoh', 'T_chnsarea', 'nsegs', 'f_selected');
    end
    load(ciCohPhasefile, 'deltaphis_allChnsSegs', 'ciCoh', 'T_chnsarea', 'nsegs', 'f_selected', 'psedociCohs');
    if ~exist('psedociCohs','var') || size(psedociCohs, 4) < shuffleN_psedoTest
        files = dir(fullfile(inputfolder_Rest, ['*_' pdcond '_*.mat']));
        [lfpsegs, fs, T_chnsarea]= seg2ShortSegments(files, twin);
        % remove unused chns
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(~removedChns_mask, :, :);
        T_chnsarea = T_chnsarea(~removedChns_mask, :);
        clear files removedChns_mask
        

        psedociCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile);
        load(ciCohPhasefile, 'psedociCohs');
        clear t_minmax_reach lfpsegs
    end
    nshuffle = size(psedociCohs, 4);
    
    % extract sigciCoh
    [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh);
    
    
    % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
    [sigciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCoh, deltaphis_allChnsSegs, T_chnsarea);
    [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, sigciCoh_flatten, deltaphis_flatten, removed_chns);
    
    
    % plot and save ciCoh Histogram image
    titlename = [animal ' Rest-'  pdcond ', nsegs = ' num2str(nsegs) ', nshuffle= ' num2str(nshuffle)];
    plot_ciCohHistogram(ciCoh_flatten_used, chnPairNames_used, f_selected, titlename);
    saveimgname = [animal '_' pdcond '.' image_type];
    saveas(gcf, fullfile(savefolder, saveimgname), image_type);
    clear titlename  saveimgname
    
    
    % rose histogram of deltaphis_allChnsTrials
    if runRose
        titlename_prefix = [animal ' Rest-'  pdcond];
        subtitlename = '';
        savefile_prefix = [animal 'trialPhaseDiff'];
        savefile_suffix = '';
        plotsave_deltaphirose(deltaphis_flatten_used, ciCoh_flatten_used, chnPairNames_used, f_selected, titlename_prefix, subtitlename, subpdsavefolder, savefile_prefix, savefile_suffix, image_type);
        clear titlename_prefix subtitlename savefile_prefix savefile_suffix
    end
    
    
    clear pdcond
    clear ciCohPhasefile
end




function [lfpdata, fs, T_chnsarea]= seg2ShortSegments(files, twin)
if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpdata = [];
    fs = [];
    T_chnsarea = [];
    
    return;
end

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