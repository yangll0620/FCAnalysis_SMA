function m4_fs500Hz_uNHP_imCohChanges_compCond(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_fs500Hz_uNHP_imCohChanges_compCond('Kitty', 'compCi_str', 1, 'ei_str', 1)
%           m4_fs500Hz_uNHP_imCohChanges_compCond('Kitty', 'compCi_str', 1, 'ei_str', 1, 'shuffleN_psedoTest', 10, 'newRun', true)
%   
%   Input:
%       animal
%
%       Name-Value:
%           compCi_str - compare cond start index
%           compCi_end - compare cond end index
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
addParameter(p, 'compCi_str', 1, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'compCi_end', [], @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', [], @isscalar);
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
compCi_str = p.Results.compCi_str;
compCi_end = p.Results.compCi_end;
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
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
EventPhases = SKT_eventPhases_extract(animal, 'codesavefolder', savecodefolder);
if isempty(ei_end)
    ei_end = length(EventPhases);
end
[tbl_compConds, ~]= compCondEvents_extract(animal);
if isempty(compCi_end)
    compCi_end = height(tbl_compConds);
end
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

for compCi = compCi_str: compCi_end
    basepd = tbl_compConds.baseCond{compCi};
    comppd = tbl_compConds.compCond{compCi};
    
    
    % extract lfptrials
    files_base = dir(fullfile(inputfolder_SKTlfp, ['*_' basepd '_*.mat']));
    if isempty(files_base)
        clear basepd comppd align2 t_AOI align2name
        clear files_base
        continue;
    end
    files_comp = dir(fullfile(inputfolder_SKTlfp, ['*_' comppd '_*.mat']));
    if isempty(files_comp)
        clear basepd comppd align2 t_AOI align2name
        clear files_base files_comp
        continue;
    end  

    for ei = ei_str : ei_end
        event = EventPhases{ei};
        [align2_base, t_AOI, align2name_base] = SKT_EventPhase_align2_tAOI_extract(event, animal, basepd, 'codesavefolder', savecodefolder);
        [align2, ~, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, comppd, 'codesavefolder', savecodefolder);
        
        ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix  '_b' basepd '--' comppd '_' event '_align2' align2name '.mat']);

        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ciCohChangesfile, 'file') || newRun)
            load(fullfile(inputfolder_SKTciCoh, [inputfile_cicoh_prefix basepd '_' event '_align2' align2name_base '.mat']), 'ciCoh', 'ntrials', 'f_selected', 'T_chnsarea');
            ciCoh_base = ciCoh; clear ciCoh
            load(fullfile(inputfolder_SKTciCoh, [inputfile_cicoh_prefix comppd '_' event '_align2' align2name '.mat']), 'ciCoh');
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
            if strcmpi(animal, 'Kitty')
                [lfptrials_base, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files_base, align2_base,t_AOI, 'codesavefolder', savecodefolder);
                [lfptrials_comp, ~, ~] = lfpseg_selectedTrials_align2(files_comp, align2,t_AOI, 'codesavefolder', savecodefolder);
            end
            if strcmpi(animal, 'Jo')
                [lfptrials_base, fs, T_chnsarea] = lfptrials_goodTrials_align2(files_base, align2_base,t_AOI, 'codesavefolder', savecodefolder);
                [lfptrials_comp, ~, ~] = lfptrials_goodTrials_align2(files_comp, align2,t_AOI, 'codesavefolder', savecodefolder);
            end
            
            
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);
            
            
            lfptrials_base = lfptrials_base(mask_chnOfI, :, :);
            lfptrials_comp = lfptrials_comp(mask_chnOfI, :, :);
            
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials_comp, lfptrials_base, fs, f_AOI, ciCohChangesfile, 'codesavefolder', savecodefolder);
            
            clear lfptrials_base fs T_chnsarea mask_chnOfI lfptrials_base lfptrials_comp
        end
        clear psedoVar

        %%% -- plot and save section --- %%%
        load(ciCohChangesfile, 'ciCohChanges', 'psedoiCohChanges', 'f_selected', 'ntrials', 'T_chnsarea')
        
        % extract sigciCoh
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        
        % plot and save ciCoh Histogram image
        nshuffle = size(psedoiCohChanges,4);
        titlename = [imgtitle_prefix ':' comppd '/' basepd '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name];
        plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, titlename, 'histClim', [-1 1],...
            'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1]);
        saveimgname = [saveimg_prefix '-' event '-b' basepd  '-' comppd '_align2' char(align2) '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        close(gcf)
        
        %%% clear
        clear event ciCohChangesfile
        clear ciCohChanges psedoiCohChanges f_selected ntrials
        clear sigciCohChanges sigciCohChanges_flatten chnPairNames
        clear titlename saveimgname nshuffle
        clear align2 t_AOI align2name align2name_base align2_base
    end
    
    clear basepd comppd 
    clear files_base files_comp
    
end


