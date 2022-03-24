function m4_fs1000Hz_FreezeSKT_imCohUsingFFT_EqualDurSegnum(varargin)
% 
%   Input
%       Name-Value: 
%           'matchSKT' - tag for matchSKT, true or false (default)
%           'ntrialsUsed' - ntrials used for cicoh calculation, default = 100(matchSKT should be false)

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
addpath(genpath(fullfile(codefolder,'connAnalyTool')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);


% parse params
p = inputParser;
addParameter(p, 'matchSKT', false, @(x)isscalar(x)&&islogical(x));
addParameter(p, 'ntrialsUsed', 100, @(x)isscalar(x)&&isnumeric(x));
parse(p,varargin{:});
matchSKT = p.Results.matchSKT;
ntrialsUsed = p.Results.ntrialsUsed;


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


ciCohPhasefile_prefix =[animal '_Freeze_ciCoh'];
if ~matchSKT
   ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_ntrial' num2str(ntrialsUsed)];
else
    ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_matchSKT'];
end

%%  input setup
inputfolder_Freeze = fullfile(codecorresParentfolder, 'm3_fs1000Hz_freezeSKTData_EpisodeExtract');
inputfolder_SKT = fullfile(codecorresParentfolder,'..' ,'m4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');


twin = 0.2;

image_type = 'tif';

f_AOI = [8 40];

pdcond = 'moderate';

shuffleN_psedoTest = 500;


%% Code start here

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
notAOI_chns = notInterested_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns

    
subpdsavefolder = fullfile(savefolder, pdcond);
if ~exist(subpdsavefolder, 'dir')
    mkdir(subpdsavefolder);
end


% load(and extract) ciCohPhasefile
ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond  '.mat']);
if(~exist(ciCohPhasefile, 'file'))

    files = dir(fullfile(inputfolder_Freeze, ['*_' pdcond '_*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= seg2ShortSegments(files, twin);
    
    save(ciCohPhasefile, 'T_chnsarea', 'combFreeTypes');
    
    
    ciCohs = struct();
    
    % for each freeze type calculate cicoh and psedoCicoh
    for frTi = 1 : length(combFreeTypes)
        freezType = combFreeTypes{frTi};
        disp(['freezType = ' freezType])
        
        lfpsegs = lfpsegs_freeze.(freezType);
        
        % remove unused chns
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(~removedChns_mask, :, :);
        T_chnsarea = T_chnsarea(~removedChns_mask, :);
        clear files removedChns_mask
        
        if matchSKT
            % match the trial number in rest and SKT
            load(fullfile(inputfolder_SKT, [animal ' ciCohPhasefile_' pdcond '_earlyReach_align2ReachOnset.mat']), 'ntrials')
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrials);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
            clear ntrials nsegs randomSKTInds
        else
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrialsUsed);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
        end
        
        
        %  extract and save deltaphis_allChnsTrials and cicoh
        [~, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
        
        
        % save
        nsegs = size(lfpsegs, 3);
        ciCohs.(freezType) = ciCoh;
        save(ciCohPhasefile, 'ciCohs', 'nsegs', 'f_selected', '-append');
        
        
        %  extract and save psedociCohs
        psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, freezType, 'codesavefolder', savecodefolder);
        
        clear lfpsegs nsegs ciCoh f_selected freezType
    end
    
    clear files lfpsegs_freeze fs T_chnsarea combFreeTypes cicohs
end


load(ciCohPhasefile,'combFreeTypes', 'psedociCohs');
if ~exist('psedociCohs','var')
    files = dir(fullfile(inputfolder_Freeze, ['*_' pdcond '_*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea]= seg2ShortSegments(files, twin);
    
    % for each freeze type calculate cicoh and psedoCicoh
    for frTi = 1 : length(combFreeTypes)
        freezType = combFreeTypes{frTi};
        disp(['freezType = ' freezType])
        
        lfpsegs = lfpsegs_freeze.(freezType);
        
        % remove unused chns
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(~removedChns_mask, :, :);
        T_chnsarea = T_chnsarea(~removedChns_mask, :);
        clear files removedChns_mask
        
        if matchSKT
            % match the trial number in rest and SKT
            load(fullfile(inputfolder_SKT, [animal ' ciCohPhasefile_' pdcond '_earlyReach_align2ReachOnset.mat']), 'ntrials')
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrials);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
            clear ntrials nsegs randomSKTInds
        else
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrialsUsed);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
        end
        
        %  extract and save psedociCohs
        psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, freezType);
        
        clear lfpsegs nsegs ciCoh f_selected freezType
    end
    
    clear files lfpsegs_freeze fs 
end

load(ciCohPhasefile,  'combFreeTypes', 'psedociCohs');
for frTi = 1 : length(combFreeTypes)
    freezType = combFreeTypes{frTi};
    disp(['freezType = ' freezType])
    
    if ~isfield(psedociCohs, freezType) || size(psedociCohs.(freezType), 4) < shuffleN_psedoTest
        files = dir(fullfile(inputfolder_Freeze, ['*_' pdcond '_*.mat']));
        [lfpsegs_freeze, fs, T_chnsarea]= seg2ShortSegments(files, twin);
        lfpsegs = lfpsegs_freeze.(freezType);
        
        % remove unused chns
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(~removedChns_mask, :, :);
        T_chnsarea = T_chnsarea(~removedChns_mask, :);
        clear files removedChns_mask
        
        if matchSKT
            % match the trial number in rest and SKT
            load(fullfile(inputfolder_SKT, [animal ' ciCohPhasefile_' pdcond '_earlyReach_align2ReachOnset.mat']), 'ntrials')
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrials);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
            clear ntrials nsegs randomSKTInds
        else
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrialsUsed);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
            clear randomSKTInds
        end
        
        %  extract and save psedociCohs
        psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, freezType);
        
        clear lfpsegs_freeze  lfpsegs
    end
end


% plot and save
load(ciCohPhasefile, 'ciCohs', 'T_chnsarea', 'nsegs', 'f_selected', 'combFreeTypes', 'psedociCohs');
for frTi = 1 : length(combFreeTypes)
    freezType = combFreeTypes{frTi};
    [sigciCoh]= sigciCoh_extract(psedociCohs.(freezType), ciCohs.(freezType));
    
    [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);
    
    nshuffle = size(psedociCohs.(freezType), 4);
    titlename = [animal ' Freeze -'  freezType ', nsegs = ' num2str(nsegs) ', nshuffle= ' num2str(nshuffle)];
    plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, titlename);
    saveimgname = [animal '_Freeze' freezType '.' image_type];
    saveas(gcf, fullfile(savefolder, saveimgname), image_type);
    
    clear titlename  saveimgname 
    clear freeType sigciCoh sigciCoh_flatten chnPairNames
end




   
function [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= seg2ShortSegments(files, twin)
if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpsegs_freeze = [];
    fs = [];
    T_chnsarea = [];
    
    return;
end

t_ThreFreeze = 3;
optFreezeTypes = optFreezeTypes_extract();

combFreeTypes = {'InitFreeze', 'ReachFreeze', 'ManipuFreeze'}; % combined {'freeze during React-Reach'}  and  {'freeze during Reach'} 

lfpsegs_freeze = struct();
for frTi = 1 : length(combFreeTypes)
    lfpsegs_freeze.(combFreeTypes{frTi}) = [];
end

for fi = 1: length(files)
    loadfilename = files(fi).name;
    load(fullfile(files(fi).folder, loadfilename), 'lfpdata', 'fs_lfp', 'T_chnsarea', 'freezStruct', 'selectedTrials');
    
    % check fs and T_chnsarea same across files
    if ~exist('fs_unit', 'var')
        fs_unit = fs_lfp;
    else
        if(fs_unit ~= fs_lfp)
            disp(['fs ~= fs_lfp for file ' loadfilename])
            continue;
        end
    end
    
    if ~exist('T_chnsarea_unit', 'var')
        T_chnsarea_unit = T_chnsarea;
    else
        if(~isequal(T_chnsarea_unit, T_chnsarea))
            disp(['T_chnsarea_unit ~= T_chnsarea for file ' loadfilename])
        end
    end
    
    nwin = round(twin * fs_unit);
    
    freezEpisodes = freezStruct.freezEpisodes;
    for frzi = 1 : length(freezEpisodes) 
        
        %%%  extract freeze lfp data: seglfp
        tri = freezEpisodes{frzi}.triali;
        if ~selectedTrials(tri)
            clear tri
            continue;
        end
        
        t_str = freezEpisodes{frzi}.freezeTPhaseS(1);
        t_end = freezEpisodes{frzi}.freezeTPhaseS(2);
        if t_end - t_str < t_ThreFreeze
            clear tri t_str t_end
            continue;
        end
        idx_str = round(t_str * fs_lfp);
        idx_end = round(t_end * fs_lfp);
        seglfp  = lfpdata{tri}(:, idx_str: idx_end);

        
        
        %%% segment freeze lfpdata into short leg: shortlfp
        shortlfp = [];
        len = size(seglfp, 2);
        shortSegn = floor(len / nwin);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nwin + 1;
            endi = shortsegi * nwin;
            shortlfp = cat(3, shortlfp, seglfp(:, stri : endi));
            clear stri endi
        end
        
        freezeType = freezEpisodes{frzi}.freezeType;
        frTi = find(strcmp(freezeType, optFreezeTypes));
        if frTi ==3 || frTi == 4
            frTi = frTi -1;
        end
        lfpsegs_freeze.(combFreeTypes{frTi}) = cat(3, lfpsegs_freeze.(combFreeTypes{frTi}), shortlfp);
        
        clear tri t_str t_end freezeType idxFreeT
        clear len shortSegn shortsegi
    end   
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;