function m4_fs500Hz_FreezeSKT_earlyMiddLateFreeze_imCohUsingFFT(freezePhase, varargin)
% 
%   Usage
%       m4_fs500Hz_FreezeSKT_earlyMiddLateFreeze_imCohUsingFFT('earlyFreeze', 'ntrialsUsed', 84);
%       m4_fs500Hz_FreezeSKT_earlyMiddLateFreeze_imCohUsingFFT('middleFreeze', 'ntrialsUsed', 84);
%       m4_fs500Hz_FreezeSKT_earlyMiddLateFreeze_imCohUsingFFT('lateFreeze', 'ntrialsUsed', 84);
%
%   Input
%      freezePhase : one of {'earlyFreeze', 'middleFreeze', 'lateFreeze'}
%
%       Name-Value: 
%           'matchSKT' - tag for matchSKT, true or false (default)
%           'ntrialsUsed' - ntrials used for cicoh calculation, default = 100(matchSKT should be false)
%           'nRandom' - randomly select nRandom times, default 10

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'toolbox')));
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
addParameter(p, 'ntrialsUsed', 200, @(x)isscalar(x)&&isnumeric(x));
addParameter(p, 'shuffleN_psedoTest', 500, @(x)isscalar(x)&&isnumeric(x));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'nRandom', 10, @(x)isscalar(x)&&isnumeric(x));


parse(p,varargin{:});
matchSKT = p.Results.matchSKT;
ntrialsUsed = p.Results.ntrialsUsed;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
nRandom = p.Results.nRandom;


switch freezePhase
    case 'earlyFreeze'
        seg_tseg = 0.2;
        seg_align2 = 'startFreeze';
        seg_tdur = [0 1];
    case 'middleFreeze'
        seg_tseg = 0.2;
        seg_align2 = 'startFreeze';
        seg_tdur = [2.5 3.5];
    case 'lateFreeze'
        seg_tseg = 0.2;
        seg_align2 = 'endFreeze';
        seg_tdur = [-1 0];
    otherwise 
        disp('Wrong freeze phase');
        return
end

 


%%  input setup
inputfolder_Freeze = fullfile(codecorresParentfolder, 'm3_fs500Hz_freezeSKTData_EpisodeExtract');
if matchSKT
    inputfolder_SKT = fullfile(codecorresParentfolder, 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');
end


tseg = 0.2;

image_type = 'tif';

f_AOI = [8 40];

pdcond = 'moderate';

ignoreReachFreeze = false;
combiFreeName = 'combinedfreezTypes';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

file_prefix = [animal '-Freeze' '-' freezePhase];
ciCohPhasefile_prefix =[file_prefix '_ciCoh'];
saveimgfile_prefix = [file_prefix '_ciCohHist'];

if ~matchSKT
    if isempty(ntrialsUsed)
        ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_alltrials'];
        saveimgfile_prefix = [saveimgfile_prefix '_alltrials'];
    else
        ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_ntrial' num2str(ntrialsUsed)];
        saveimgfile_prefix = [saveimgfile_prefix '_ntrial' num2str(ntrialsUsed)];
    end
else
    ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_matchSKT'];
    saveimgfile_prefix = [saveimgfile_prefix '_matchSKT'];
end


%% Code start here

unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
notAOI_chns = notInterested_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns

% load(and extract) ciCohPhasefile
ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond  '.mat']);

%%% ----------- case no ciCohPhasefile or new Run -------- %%%
if(~exist(ciCohPhasefile, 'file') || newRun)
    
    files = dir(fullfile(inputfolder_Freeze, ['*' pdcond '*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= seg2ShortSegments_wPhase(files, 'tseg', seg_tseg, 'align2', seg_align2, 'tdur', seg_tdur);
    
    save(ciCohPhasefile, 'combFreeTypes');
    clear files
    
    % combined lfpsegs
    fnames = fields(lfpsegs_freeze);
    lfpsegs = [];
    for fname = fnames'
        lfpsegs = cat(3, lfpsegs, lfpsegs_freeze.(fname{1}));
    end
    clear fnames fname lfpsegs_freeze
    
    % remove unused chns
    removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
    lfpsegs = lfpsegs(~removedChns_mask, :, :);
    T_chnsarea = T_chnsarea(~removedChns_mask, :);
    clear files removedChns_mask
    
    
    ciCoh_ntimes = [];
    for nRandi = 1 : nRandom
        
        % select matchSKT or ntrialsUsed trials
        if matchSKT
            % match the trial number in rest and SKT
            load(fullfile(inputfolder_SKT, [animal ' ciCohPhasefile_' pdcond '_earlyReach_align2ReachOnset.mat']), 'ntrials')
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrials);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
            clear ntrials nsegs randomSKTInds
        else
            if ~isempty(ntrialsUsed)
                nsegs = size(lfpsegs, 3);
                randomSKTInds =  randsample(nsegs,ntrialsUsed);
                lfpsegs = lfpsegs(:, :, randomSKTInds);
                clear nsegs randomSKTInds
            end
        end
        
        
        %  extract and save deltaphis_allChnsTrials and cicoh
        [~, ciCoh_1time, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
        
        ciCoh_ntimes = cat(4, ciCoh_ntimes, ciCoh_1time);
        
        nseg = size(lfpsegs, 3);
        
        clear ciCoh_1time
    end
    ciCoh = squeeze(mean(ciCoh_ntimes, 4));
    
    % save
    nsegs.(combiFreeName) = nseg;
    ciCohs.(combiFreeName) = ciCoh;
    save(ciCohPhasefile, 'ciCohs', 'nsegs', 'f_selected', 'T_chnsarea', '-append');
    
    
    clear lfpsegs ciCoh f_selected nsegs T_chnsarea
end



%%%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
psedoVar = whos('-file',ciCohPhasefile, 'psedociCohs');
if isempty(psedoVar) 

    files = dir(fullfile(inputfolder_Freeze, ['*' pdcond '*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea, ~]= seg2ShortSegments_wPhase(files, 'tseg', seg_tseg, 'align2', seg_align2, 'tdur', seg_tdur);
    clear files
    
    % combined lfpsegs
    fnames = fields(lfpsegs_freeze);
    lfpsegs = [];
    for fname = fnames'
        lfpsegs = cat(3, lfpsegs, lfpsegs_freeze.(fname{1}));
    end
    clear fnames fname lfpsegs_freeze
    
    % remove unused chns
    removedChns_mask = cellfun(@(x) contains(x, removed_chns), T_chnsarea.brainarea);
    lfpsegs = lfpsegs(~removedChns_mask, :, :);
    T_chnsarea = T_chnsarea(~removedChns_mask, :);
    clear files removedChns_mask
    
    
    
    % select matchSKT or ntrialsUsed trials
    if matchSKT
        % match the trial number in rest and SKT
        load(fullfile(inputfolder_SKT, [animal ' ciCohPhasefile_' pdcond '_earlyReach_align2ReachOnset.mat']), 'ntrials')
        nsegs = size(lfpsegs, 3);
        randomSKTInds =  randsample(nsegs,ntrials);
        lfpsegs = lfpsegs(:, :, randomSKTInds);
        clear ntrials nsegs randomSKTInds
    else
        if ~isempty(ntrialsUsed)
            nsegs = size(lfpsegs, 3);
            randomSKTInds =  randsample(nsegs,ntrialsUsed);
            lfpsegs = lfpsegs(:, :, randomSKTInds);
            clear nsegs randomSKTInds
        end
    end
    
    psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, combiFreeName, 'codesavefolder', savecodefolder);
    
    clear lfpsegs fs f_AOI
end

%% plot and save
load(ciCohPhasefile, 'ciCohs', 'T_chnsarea', 'nsegs', 'f_selected', 'psedociCohs');

[sigciCoh]= sigciCoh_extract(psedociCohs.(combiFreeName), ciCohs.(combiFreeName));
[sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);

nshuffle = size(psedociCohs.(combiFreeName), 4);
imgtitle_prefix = [file_prefix];
saveimgfile_prefix = [saveimgfile_prefix];

show_freqNumLabel = false;
if strcmp(freezePhase, 'lateFreeze')
    show_freqNumLabel = true;
end

titlename = [imgtitle_prefix  ', nsegs = ' num2str(nsegs.(combiFreeName)) ', nshuffle= ' num2str(nshuffle)];
plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, titlename, ...
    'fig_width', 400, 'fig_height', 200, 'margin_inner', [5 5 50 40],...
    'codesavefolder', savecodefolder, ...
    'show_titlename', false, 'show_xlabel', show_freqNumLabel, 'show_yticklabels', false);
saveimgname = [saveimgfile_prefix '.' image_type];
saveas(gcf, fullfile(savefolder, saveimgname), image_type);

close all
end


   
function [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= seg2ShortSegments_wPhase(files, varargin)
%
%
%
%   Input:
%       files: 
%    
%       Name-Value: 
%           
%           tseg: the seg duration in second, default = 0.2
%
%           align2: the time point 0, default '' means no align2, {'startFreeze', 'endFreeze'}
%   
%           tdur: the time duration respect to align2, only needed when align2 not ''
%
%           tstart: start used time point, default 0
%
if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpsegs_freeze = [];
    fs = [];
    T_chnsarea = [];
    combFreeTypes = [];
    
    return;
end

% parse params
p = inputParser;
addParameter(p, 'tseg', 0.2, @isscalar);
addParameter(p, 'align2', '', @isstr);
addParameter(p, 'tdur', [], @isvector);
addParameter(p, 'tstart', 0, @isvector);
parse(p,varargin{:});
tseg = p.Results.tseg;
align2 = p.Results.align2;
tdur = p.Results.tdur;
tstart = p.Results.tstart;
if isempty(align2)
   tdur = []; 
end



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
    
    nseg = round(tseg * fs_unit);
    
    freezEpisodes = freezStruct.freezEpisodes;
    for frzi = 1 : length(freezEpisodes) 
        
        %%%  extract freeze lfp data: seglfp
        tri = freezEpisodes{frzi}.triali;
        if ~selectedTrials(tri)
            clear tri
            continue;
        end
        if isempty(tdur)
            t_str = freezEpisodes{frzi}.freezeTPhaseS(1) + tstart;
            t_end = freezEpisodes{frzi}.freezeTPhaseS(2);
        else
            switch align2
                case 'startFreeze'
                    t0 = freezEpisodes{frzi}.freezeTPhaseS(1);
                case 'endFreeze'
                    t0 = freezEpisodes{frzi}.freezeTPhaseS(2);
            end
            t_str = tdur(1) + t0;
            t_end = tdur(2) + t0;
        end
        
        if t_end - t_str < tseg
            clear tri t_str t_end
            continue;
        end
        
        idx_str = round(t_str * fs_lfp);
        idx_end = round(t_end * fs_lfp);
        seglfp  = lfpdata{tri}(:, idx_str: idx_end);

        
        %%% segment freeze lfpdata into short leg: shortlfp
        shortlfp = [];
        len = size(seglfp, 2);
        shortSegn = floor(len / nseg);
        for shortsegi = 1 : shortSegn
            stri = (shortsegi - 1)* nseg + 1;
            endi = shortsegi * nseg;
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
    
    clear nseg
end
fs = fs_unit;
T_chnsarea = T_chnsarea_unit;
end