function m4_fs500Hz_uNHP_imCohChanges_compEvent(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_fs500Hz_uNHP_imCohChanges_compEvent('Kitty', 'compEi_str', 1, 'ci_str', 1)
%           m4_fs500Hz_uNHP_imCohChanges_compEvent('Kitty', 'compEi_str', 1, 'ci_str', 1, 'shuffleN_psedoTest', 10, 'newRun', true)
%   
%   Input:
%       animal
%
%       Name-Value:
%           compEi_str - compare events start index
%           compEi_end - compare events end index
%           ei_str - event start index
%           ei_end - event end index
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500

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
addParameter(p, 'compEi_str', 1, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'compEi_end', [], @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'ci_str', 2, @isscalar);
addParameter(p, 'ci_end', [], @isscalar);
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
compEi_str = p.Results.compEi_str;
compEi_end = p.Results.compEi_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end

[~, ~, pipelinefolder, ~] = exp_subfolders();

%% save setup
[~,codefilename,~] = fileparts(codefilepath); 
savecodefolder = fullfile(codefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', codefilename);
[savefolder, ~] = code_corresfolder(savecodefolder, false, false);

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

ciCohChangesfile_prefix =[animal '-ciCohChanges'];
saveimg_prefix = [animal '-ciCohChanges'];
imgtitle_prefix = [animal '-FCChanges'];

%%  input setup
if strcmpi(animal, 'Kitty')
    cicohfolder = 'm3_fs500Hz_uNHP_histFTLagPhase';
    lfpfolder = 'm2_segSKTData_SelectTrials_chnOfI';
    inputfile_cicoh_prefix = [animal ' ciCohPhasefile_8-40Hz_'];
end
if strcmpi(animal, 'Jo')
    cicohfolder = 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP';
    lfpfolder = 'm2_SKTData_SelectTrials';
    inputfile_cicoh_prefix = [animal ' ciCohPhasefile_'];
end

inputfolder_SKTciCoh = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', cicohfolder);
inputfolder_SKTlfp = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', lfpfolder);

image_type = 'tif';

f_AOI = [8 40];



%% Code start here
cond_cell = cond_cell_extract(animal);
if isempty(ci_end)
    ci_end = length(cond_cell);
end
[~, tbl_compEvents]= compCondEvents_extract(animal);
if isempty(compEi_end)
    compEi_end = height(tbl_compEvents);
end
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

for compEi = compEi_str: compEi_end
    baseEvent = tbl_compEvents.baseEvent{compEi};
    compEvent = tbl_compEvents.compEvent{compEi};

    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};
        
        [align2_base, t_AOI, align2name_base] = SKT_EventPhase_align2_tAOI_extract(baseEvent, animal, pdcond, 'codesavefolder', savecodefolder);
        [align2, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(compEvent, animal, pdcond, 'codesavefolder', savecodefolder);
        
        ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix '_' pdcond '_b' baseEvent '--' compEvent  '_align2' align2name '.mat']);

        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ciCohChangesfile, 'file') || newRun)
            load(fullfile(inputfolder_SKTciCoh, [inputfile_cicoh_prefix pdcond '_' baseEvent '_align2' align2name_base '.mat']), 'ciCoh', 'ntrials', 'f_selected', 'T_chnsarea');
            ciCoh_base = ciCoh; clear ciCoh
            load(fullfile(inputfolder_SKTciCoh, [inputfile_cicoh_prefix pdcond '_' compEvent '_align2' align2name '.mat']), 'ciCoh');
            ciCoh_comp = ciCoh; clear ciCoh
            
            ciCohChanges = ciCoh_comp - ciCoh_base;
            if strcmpi(animal, 'Jo')
                mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                T_chnsarea = T_chnsarea(mask_chnOfI, :);
                ciCohChanges = ciCohChanges(mask_chnOfI, mask_chnOfI, :);
                clear mask_chnOfI
            end
            
            save(ciCohChangesfile, 'ciCohChanges', 'ntrials', 'f_selected', 'T_chnsarea')
            
            clear ciCoh_base ciCoh_comp 
            clear ciCohChanges ntrials f_selected
        end
            
        
        %%%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
        psedoVar = whos('-file',ciCohChangesfile, 'psedoiCohChanges');
        if isempty(psedoVar) || psedoVar.size(4)< shuffleN_psedoTest
            
            % extract lfptrials
            lfpfiles = dir(fullfile(inputfolder_SKTlfp, ['*_' pdcond '_*.mat']));
            if isempty(lfpfiles)
                clear lfpfiles
                continue;
            end
            
            if strcmpi(animal, 'Kitty')
                [lfptrials_base, fs, T_chnsarea] = lfpseg_selectedTrials_align2(lfpfiles, align2_base,t_AOI, 'codesavefolder', savecodefolder);
                [lfptrials_comp, ~, ~] = lfpseg_selectedTrials_align2(lfpfiles, align2,t_AOI, 'codesavefolder', savecodefolder);
            end
            if strcmpi(animal, 'Jo')
                [lfptrials_base, fs, T_chnsarea] = lfptrials_goodTrials_align2(lfpfiles, align2_base,t_AOI, 'codesavefolder', savecodefolder);
                [lfptrials_comp, ~, ~] = lfptrials_goodTrials_align2(lfpfiles, align2,t_AOI, 'codesavefolder', savecodefolder);
            end
            
            
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);
            
            
            lfptrials_base = lfptrials_base(mask_chnOfI, :, :);
            lfptrials_comp = lfptrials_comp(mask_chnOfI, :, :);
            
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials_comp, lfptrials_base, fs, f_AOI, ciCohChangesfile, 'codesavefolder', savecodefolder);
            
            clear lfpfiles lfptrials_base fs T_chnsarea mask_chnOfI lfptrials_base lfptrials_comp
        end
        clear psedoVar

        %%% -- plot and save section --- %%%
        load(ciCohChangesfile, 'ciCohChanges', 'psedoiCohChanges', 'f_selected', 'ntrials', 'T_chnsarea')
        
        % extract sigciCoh
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        
        % plot and save ciCoh Histogram image
        titlename = [imgtitle_prefix ':' pdcond '-' compEvent '/' baseEvent];
        plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, titlename, 'histClim', [-1 1],...
            'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1]);
        saveimgname = [saveimg_prefix '-' pdcond '-b' baseEvent  '-' compEvent '_align2' char(align2) '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        close(gcf)
        
        %%% clear
        clear pdcond ciCohChangesfile
        clear ciCohChanges psedoiCohChanges f_selected ntrials
        clear sigciCohChanges sigciCohChanges_flatten chnPairNames
        clear titlename saveimgname nshuffle
        clear align2 t_AOI align2name align2name_base align2_base
    end
    
    clear compEvent baseEvent     
end


