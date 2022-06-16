function m7_imCohChanges_chnOfI_unifiedNHP(animal, varargin)
% plot cicoh Histogram, frequency time lag and phase of interested channels
%   Input:
%       Name-Value: 
%           animal
%           ei_str - event start index
%           ci_str - condition start index

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


% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
SKTSubfolder = 'SKT';
if strcmpi(animal, 'Kitty')
    SKTSubfolder = 'SKT_SegV';
end
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);

%% Input setup
inputfolder = fullfile(codecorresParentfolder, 'm6_imCohChangesUsingFFT_basedNormal_unifiedNHP');
image_type = 'tif';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code start here

cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();
chnsOfI = chnOfInterest_extract(animal, 'codesavefolder', savecodefolder);

for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    for ei = 1 : length(EventPhases)
        ephase = EventPhases{ei};
        [~, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(ephase, animal, pdcond);
        filename = [animal ' ciCohChangesfile_' pdcond '_' ephase '_align2' align2name '.mat'];
        if ~exist(fullfile(inputfolder, filename), 'file')
            continue;
        end
        
        load(fullfile(inputfolder, filename), 'f_selected', 'ntrials','psedoiCohChanges', 'ciCohChanges', 'T_chnsarea');
        
        
        % select the data of chnOfI
        mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
        T_chnsarea = T_chnsarea(mask_chnOfI, :);
        ciCohChanges = ciCohChanges(mask_chnOfI, mask_chnOfI, :);
        psedoiCohChanges = psedoiCohChanges(mask_chnOfI, mask_chnOfI, :, :);
        clear mask_chnOfI
        
        
        % extract sigciCoh
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        
        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        
        
        % plot and save ciCohChanges Histogram image
        nshuffle = size(psedoiCohChanges,4);
        titlename = [animal ' FCChanges '  pdcond '-'  ephase '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
        plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, round(f_selected), titlename, [-1 1], ...
            'fig_width', 1000, 'fig_height', 250);
        saveimgname = [animal 'FCChangesChoI_' ephase '_' pdcond '_align2' align2name '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        clear titlename  saveimgname nshuffle
        close all

        clear ephase align2name filename t_AOI
        clear('f_selected', 'ntrials','psedoiCohChanges', 'ciCohChanges');
    end
end
