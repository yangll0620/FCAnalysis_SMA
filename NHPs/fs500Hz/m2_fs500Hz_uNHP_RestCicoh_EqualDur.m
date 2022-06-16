function m2_fs500Hz_uNHP_RestCicoh_EqualDur(animal, varargin)
% plot cicoh Histogram, frequency time lag and phase of interested channels
%
%   Example Usage:
%           m2_fs500Hz_uNHP_RestCicoh_EqualDur('Kitty', 'matchSKT', true, 'ci_str', 1)
%   
%   Input:
%       animal
%
%       Name-Value:
%           matchSKT - tag for match SKT or not (default true)l;
%           nsegsUsed - number of segs used for calc rest cicoh (default 200 trials, used only when matchSKT = false)
%           ci_str - event start index
%           ci_end - condition end index
%           newRun - true or false(default), for running new or not

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


% parse params
p = inputParser;
addParameter(p, 'ci_str', 1, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'ci_end', [], @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'nsegsUsed', 200, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'matchSKT', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end

parse(p,varargin{:});
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
matchSKT = p.Results.matchSKT;
if ~matchSKT
    nsegsUsed = p.Results.nsegsUsed;
end
disp('p.Results =  ' )
p.Results


[~, ~, pipelinefolder, ~] = exp_subfolders();


%% save setup
[~,codefilename,~] = fileparts(codefilepath);

savecodefolder = fullfile(codefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', codefilename);
[savefolder, ~] = code_corresfolder(savecodefolder, false, false);

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

ciCohPhasefile_prefix =[animal '_Rest_ciCohPhasefile'];
if matchSKT
    ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_matchSKT'];
    saveimg_prefix = [animal '_Rest' '_matchSKT'];
else
    ciCohPhasefile_prefix = [ciCohPhasefile_prefix '_nSegs' num2str(nsegsUsed)];
    saveimg_prefix = [animal '_Rest' '_nSegs' num2str(nsegsUsed)];
end

%%  input setup
switch animal
    case 'Kitty'
        restdatafolder = 'm1_restData_rmChns_avgArea';
        if matchSKT
            SKTdatafolder = 'm3_fs500Hz_uNHP_histFTLagPhase';
        end
    case 'Jo'
        restdatafolder = 'm3_restData_rmChns_avgArea';
        if matchSKT
            SKTdatafolder = 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP';
        end
end
inputcodefolder = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'Rest', restdatafolder);
[inputfolder, ~] = code_corresfolder(inputcodefolder, false, false);

if matchSKT
    inputfolder_SKT = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', SKTdatafolder);
end

twin = 0.2;

image_type = 'tif';

f_AOI = [8 40];


%% Code start here
disp('savefolder: ')
disp(savefolder);

cond_cell = cond_cell_extract(animal);
if isempty(ci_end)
    ci_end = length(cond_cell);
end
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
for ci = ci_str : ci_end
    pdcond = cond_cell{ci};
    
    disp([animal '-' pdcond])
    
    ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' pdcond  '.mat']);
    
    if matchSKT
        fileSKT = dir(fullfile(inputfolder_SKT, [animal ' *ciCohPhasefile*' pdcond '*earlyReach_align2ReachOnset.mat']));
        if length(fileSKT) ~= 1
            disp('fileSKT not only one')
            disp(inputfolder_SKT);
            continue;
        end
        load(fullfile(inputfolder_SKT, fileSKT.name), 'ntrials');
        nsegsUsed = ntrials;
        clear ntrials
    end
    
    %%% ----------- case no ciCohPhasefile or new Run -------- %%%
    if(~exist(ciCohPhasefile, 'file') || newRun)
        files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
        if isempty(files)
            clear files
            continue;
        end
        
        % lfpseg
        [lfpsegs, fs, T_chnsarea]= seg2ShortSegments(files, twin);
        
        % reandomly selected nsegsUsed
        nsegs = size(lfpsegs, 3);
        randomInds =  randsample(nsegs,nsegsUsed);
        lfpsegs = lfpsegs(:, :, randomInds);
        clear randomInds
        
        % extract data of chns of AOI
        mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(mask_chnOfI, :, :);
        T_chnsarea = T_chnsarea(mask_chnOfI, :);
        
        
        %  extract and save deltaphis_allChnsTrials and cicoh
        [deltaphis_allChnsSegs, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI, 'codesavefolder', savecodefolder);
        nsegs = size(lfpsegs, 3);
        save(ciCohPhasefile, 'deltaphis_allChnsSegs', 'ciCoh', 'T_chnsarea', 'nsegs', 'f_selected');
        
        clear files lfpsegs fs T_chnsarea mask_chnOfI
        clear deltaphis_allChnsSegs ciCoh f_selected nsegs
    end
       
    
    %%% ----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
    pseciCohsVar = whos('-file',ciCohPhasefile, 'psedociCohs');
    if isempty(pseciCohsVar) || pseciCohsVar.size(4)< shuffleN_psedoTest
        files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
        if isempty(files)
            clear files
            continue;
        end
        
        % lfpseg
        [lfpsegs, fs, T_chnsarea]= seg2ShortSegments(files, twin);
        
        % reandomly selected nsegsUsed
        nsegs = size(lfpsegs, 3);
        randomInds =  randsample(nsegs,nsegsUsed);
        lfpsegs = lfpsegs(:, :, randomInds);
        clear randomInds nsegs
        
        % extract data of chns of AOI
        mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
        lfpsegs = lfpsegs(mask_chnOfI, :, :);
        T_chnsarea = T_chnsarea(mask_chnOfI, :);
        
        %  extract and save psedociCohs
        psedociCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCohPhasefile);
        
        clear files lfpsegs fs T_chnsarea mask_chnOfI
        clear deltaphis_allChnsTrials ciCoh f_selected nsegs
    end
    clear pseciCohsVar
    
    
    %%% -- plot section --- %%%
    load(ciCohPhasefile,  'ciCoh', 'T_chnsarea', 'nsegs', 'f_selected', 'psedociCohs');
    
    % extract sigciCoh
    [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh);
    
    
    % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
    [sigciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCoh, T_chnsarea, 'codesavefolder',savecodefolder);
    
    
    % plot and save ciCoh Histogram image
    nshuffle = size(psedociCohs,4);
    titlename = [animal ' Rest-'  pdcond ', nseg = ' num2str(nsegs) ' nshuffle= ' num2str(nshuffle)];
    plot_ciCohHistogram(sigciCoh_flatten, chnPairNames, f_selected, titlename, ...
                'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder);
    saveimgname = [saveimg_prefix  '_' pdcond '.' image_type];
    saveas(gcf, fullfile(savefolder, saveimgname), image_type);
    clear titlename  saveimgname nshuffle
        
    
    clear pdcond ciCohPhasefile
    clear('ciCoh', 'T_chnsarea', 'nsegs', 'f_selected', 'psedociCohs');
    clear sigciCoh sigciCoh_flatten chnPairNames
end
close all



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