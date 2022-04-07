function m4_fs500Hz_FreezeSKT_imCohUsingFFT(varargin)
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
addParameter(p, 'ntrialsUsed', [], @(x)isscalar(x)&&isnumeric(x));
parse(p,varargin{:});
matchSKT = p.Results.matchSKT;
ntrialsUsed = p.Results.ntrialsUsed;


%%  input setup
inputfolder_Freeze = fullfile(codecorresParentfolder, 'm3_fs500Hz_freezeSKTData_EpisodeExtract');
if matchSKT
    inputfolder_SKT = fullfile(codecorresParentfolder, 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');
end


tseg = 0.2;

image_type = 'tif';

f_AOI = [8 40];

pdcond = 'moderate';

shuffleN_psedoTest = 280;
ignoreReachFreeze = true;
combiFreeName = 'combinedfreezTypes';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

file_prefix = [animal '_Freeze'];
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

    
subpdsavefolder = fullfile(savefolder, pdcond);
if ~exist(subpdsavefolder, 'dir')
    mkdir(subpdsavefolder);
end


% load(and extract) ciCohPhasefile
ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond  '.mat']);
if(~exist(ciCohPhasefile, 'file'))

    files = dir(fullfile(inputfolder_Freeze, ['*' pdcond '*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= seg2ShortSegments(files, tseg);
    
    save(ciCohPhasefile, 'T_chnsarea', 'combFreeTypes');
    
    
    ciCohs = struct();
    nsegs = struct();
    
    %%% for each freeze type calculate cicoh and psedoCicoh
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
            if ~isempty(ntrialsUsed)
                nsegs = size(lfpsegs, 3);
                randomSKTInds =  randsample(nsegs,ntrialsUsed);
                lfpsegs = lfpsegs(:, :, randomSKTInds);
                clear nsegs randomSKTInds
            end
        end
        
        
        %  extract and save deltaphis_allChnsTrials and cicoh
        [~, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
        
        
        % save
        nsegs.(freezType) = size(lfpsegs, 3);
        ciCohs.(freezType) = ciCoh;
        save(ciCohPhasefile, 'ciCohs', 'nsegs', 'f_selected', '-append');
        
        
        %  extract and save psedociCohs
        psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, freezType, 'codesavefolder', savecodefolder);
        
        clear lfpsegs ciCoh f_selected freezType
    end
    
    
    %%% calculate cicoh and psedoCicoh using combined lfp from all freeze types
    disp(['Using combined freeze type lfp'])
    lfpsegs = [];
    for frTi = 1 : length(combFreeTypes) % combine lfp
        freezType = combFreeTypes{frTi};
        if ignoreReachFreeze && strcmp(combFreeTypes{frTi}, 'ReachFreeze')
            continue;
        end
        lfpsegs = cat(3, lfpsegs, lfpsegs_freeze.(freezType));
    end
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
             
    %  extract and save deltaphis_allChnsTrials and cicoh
    [~, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
            
    % save
    nsegs.(combiFreeName) = size(lfpsegs, 3);
    ciCohs.(combiFreeName) = ciCoh;
    save(ciCohPhasefile, 'ciCohs', 'nsegs', 'f_selected', '-append');     
        
    %  extract and save psedociCohs
    psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, combiFreeName, 'codesavefolder', savecodefolder);
        
    clear lfpsegs ciCoh f_selected 
    
    %%% final clear
    clear files lfpsegs_freeze fs T_chnsarea combFreeTypes cicohs nsegs
end


load(ciCohPhasefile, 'ciCohs', 'nsegs', 'combFreeTypes', 'psedociCohs');
if ~exist('psedociCohs','var')
    files = dir(fullfile(inputfolder_Freeze, ['*_' pdcond '_*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea]= seg2ShortSegments(files, tseg);
    
    %%% for each freeze type calculate cicoh and psedoCicoh
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
    
    %%% calculate cicoh and psedoCicoh using combined lfp from all freeze types
    disp(['Using combined freeze type lfp'])
    lfpsegs = [];
    for frTi = 1 : length(combFreeTypes) % combine lfp
        freezType = combFreeTypes{frTi};
        if ignoreReachFreeze && strcmp(combFreeTypes{frTi}, 'ReachFreeze')
            continue;
        end
        lfpsegs = cat(3, lfpsegs, lfpsegs_freeze.(freezType));
    end
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
             
    %  extract and save deltaphis_allChnsTrials and cicoh
    [~, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
            
    % save
    nsegs.(combiFreeName) = size(lfpsegs, 3);
    ciCohs.(combiFreeName) = ciCoh;
    save(ciCohPhasefile, 'ciCohs', 'nsegs', 'f_selected', '-append');     
        
    %  extract and save psedociCohs
    psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, combiFreeName, 'codesavefolder', savecodefolder);
        
    clear lfpsegs ciCoh f_selected
    
    
    %%% final clear
    clear files lfpsegs_freeze fs  T_chnsarea
end

for frTi = 1 : length(combFreeTypes)
    freezType = combFreeTypes{frTi};
    disp(['freezType = ' freezType])
    
    if ~isfield(psedociCohs, freezType) || size(psedociCohs.(freezType), 4) < shuffleN_psedoTest
        files = dir(fullfile(inputfolder_Freeze, ['*_' pdcond '_*.mat']));
        [lfpsegs_freeze, fs, T_chnsarea]= seg2ShortSegments(files, tseg);
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


%%% calculate cicoh and psedoCicoh using combined lfp from all freeze types
if ~isfield(psedociCohs, combiFreeName) || size(psedociCohs.(combiFreeName), 4) < shuffleN_psedoTest
    files = dir(fullfile(inputfolder_Freeze, ['*_' pdcond '_*.mat']));
    [lfpsegs_freeze, fs, T_chnsarea]= seg2ShortSegments(files, tseg);

    disp(['Using combined freeze type lfp'])
    lfpsegs = [];
    for frTi = 1 : length(combFreeTypes) % combine lfp
        freezType = combFreeTypes{frTi};
        if ignoreReachFreeze && strcmp(combFreeTypes{frTi}, 'ReachFreeze')
            continue;
        end
        lfpsegs = cat(3, lfpsegs, lfpsegs_freeze.(freezType));
    end
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
    
    %  extract and save deltaphis_allChnsTrials and cicoh
    [~, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);
    
    % save
    nsegs.(combiFreeName) = size(lfpsegs, 3);
    ciCohs.(combiFreeName) = ciCoh;
    save(ciCohPhasefile, 'ciCohs', 'nsegs', 'f_selected', '-append');
    
    %  extract and save psedociCohs
    psedoFreeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile, combiFreeName, 'codesavefolder', savecodefolder);
    
    clear lfpsegs ciCoh f_selected
end


%% plot and save
load(ciCohPhasefile, 'ciCohs', 'T_chnsarea', 'nsegs', 'f_selected', 'combFreeTypes', 'psedociCohs');

% plot cicoh from each freeze type
for frTi = 1 : length(combFreeTypes)
    freezType = combFreeTypes{frTi};
    [sigciCoh]= sigciCoh_extract(psedociCohs.(freezType), ciCohs.(freezType));
    
    [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);
    
    nshuffle = size(psedociCohs.(freezType), 4);
    titlename = [animal ' Freeze -'  freezType ', nsegs = ' num2str(nsegs.(freezType)) ', nshuffle= ' num2str(nshuffle)];
    plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, titlename, 'fig_width', 1000, 'fig_height', 250);
    saveimgname = [saveimgfile_prefix '_Freeze' freezType '.' image_type];
    saveas(gcf, fullfile(savefolder, saveimgname), image_type);
    
    clear titlename  saveimgname 
    clear freeType sigciCoh sigciCoh_flatten chnPairNames
end

% plot cicoh from all combined freeze types
[sigciCoh]= sigciCoh_extract(psedociCohs.(combiFreeName), ciCohs.(combiFreeName));
[sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea);


nshuffle = size(psedociCohs.(combiFreeName), 4);
if ignoreReachFreeze
    imgtitle_prefix = [file_prefix '-' combiFreeName '-noReach'];
    saveimgfile_prefix = [saveimgfile_prefix '_' combiFreeName '-noReach'];
else
    imgtitle_prefix = [file_prefix '-' combiFreeName];
    saveimgfile_prefix = [saveimgfile_prefix '_' combiFreeName];
end
titlename = [imgtitle_prefix  ', nsegs = ' num2str(nsegs.(combiFreeName)) ', nshuffle= ' num2str(nshuffle)];
plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, titlename, 'fig_width', 1000, 'fig_height', 250);
saveimgname = [saveimgfile_prefix '.' image_type];
saveas(gcf, fullfile(savefolder, saveimgname), image_type);
clear titlename  saveimgname
clear freeType sigciCoh sigciCoh_flatten chnPairNames

close all
end


function plot_ciCohHistogram(ciCoh_flatten, chnPairNames, f_selected, titlename, varargin)
%
%   Inputs:
%       ciCoh_flatten:
%       chnPairNames
%       f_selected
%       titlename
%       histClim
%
%       Name-Value: 
%           'codesavefolder' - code saved folder


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [0 1], @(x) assert(isnumeric(x) && isvector(x)));
addParameter(p, 'fig_height', 600, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_width', 600, @(x) assert(isnumeric(x) && isscalar(x)));
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

fig_left = 50;
fig_bottom = 50;
fig_width = p.Results.fig_width;
fig_height = p.Results.fig_height;
histClim = p.Results.histClim;


% plot
figure;
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
imagesc(ciCoh_flatten)
colormap(jet)
set(gca, 'Position', [0.15 0.2 0.75 0.7])
[npairs, nf] = size(ciCoh_flatten);
xticks([1:nf])
xticklabels(round(f_selected))
yticks([1:npairs]);
set(gca,'YTickLabel',chnPairNames,'fontsize',12,'FontWeight','bold')
xlabel('freqs')
title(titlename, 'FontSize', 15, 'FontWeight', 'normal')
set(gca,'CLim', histClim)
colorbar



chnPair_prev = '';
for ci = 1: length(chnPairNames)
    chnPair = chnPairNames{ci};
    
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
end

   
function [lfpsegs_freeze, fs, T_chnsarea, combFreeTypes]= seg2ShortSegments(files, tseg)
if isempty(files)
    disp('files for seg2ShortSegments empty!')
    
    lfpsegs_freeze = [];
    fs = [];
    T_chnsarea = [];
    combFreeTypes = [];
    
    return;
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
        
        t_str = freezEpisodes{frzi}.freezeTPhaseS(1);
        t_end = freezEpisodes{frzi}.freezeTPhaseS(2);
        if t_end - t_str < tseg
            clear tri t_str t_end
            continue;
        end
        idx_str = round((t_str + 0.5) * fs_lfp);
        idx_end = round((t_end - 0.5) * fs_lfp);
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