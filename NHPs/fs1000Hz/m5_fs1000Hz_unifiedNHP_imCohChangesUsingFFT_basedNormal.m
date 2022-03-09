function m5_fs1000Hz_unifiedNHP_imCohChangesUsingFFT_basedNormal(animal, varargin)
%
%   Input:
%       animal
%
%       Name-Value: 
%           ei_str - event start index
%           ei_end - event end index
%           ci_str - event start index
%           ci_end - condition end index
%           runCicohHist -  true (default) or false(default)
%           runRosePlot -  true or false(default)

codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'toolbox')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));


cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();

% parse params
p = inputParser;
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', length(EventPhases), @isscalar);
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', length(cond_cell), @isscalar);
addParameter(p, 'runCicohHist', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'runRosePlot', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
disp('p.Results =  ' )
p.Results


% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
segVFolder = false;
if strcmpi(animal, 'Kitty')
    segVFolder = true;
end

if segVFolder
    SKTSubfolder = 'SKT_SegV';
    SKTlfpSubfolder = 'm2_segSKTData_SelectTrials_goodReach_fs1000Hz';
else
    SKTSubfolder = 'SKT';
    SKTlfpSubfolder = 'm2_SKTData_SelectTrials_goodReach_fs1000Hz';
end
SKTCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, 'fs1000Hz', codefilename);
[SKTCodecorresfolder, SKTCodecorresParentfolder] = code_corresfolder(SKTCodefilepath, true, false);


%% save setup
savefolder = SKTCodecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
savematfolder = fullfile(savefolder, 'mat');
saveimgfolder = fullfile(savefolder, 'images');

ciCohChangesfile_prefix =[animal '_FCChan_based'];

%%  input setup
inputfolder_ciCohSKT = fullfile(SKTCodecorresParentfolder, 'm4_fs1000Hz_unifiedNHP_chnOfI_imCohPhaseUsingFFT_EventPhase');
inputfolder_lfpSKT = fullfile(SKTCodecorresParentfolder, SKTlfpSubfolder);


image_type = 'tif';

twin = 0.2;

%% Code start here
copyfile2folder(codefilepath, savecodefolder);

if ~exist(savematfolder, 'dir')
    mkdir(savematfolder)
end
if ~exist(saveimgfolder, 'dir')
    mkdir(saveimgfolder)
end

cond_cell = cond_cell_extract(animal);

EventPhases = SKT_eventPhases_extract();
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal);

% base is normal
basepd = 'normal';
fhigh_AOIs = freqsOfInterest_extract(animal, 'ishighfreq', true, 'codesavefolder', fullfile(savefolder, 'code'));
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
for fi = 1: size(fhigh_AOIs, 1)
    f_AOI = fhigh_AOIs(fi, :);
    for ei = ei_str: ei_end
        event = EventPhases{ei};
        
        % load base ciCoh
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, basepd);
        load(fullfile(inputfolder_ciCohSKT, [animal ' ciCohPhasefile' '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz_' basepd '_' event '_align2' align2name '.mat']), 'ciCoh');
        ciCoh_base = ciCoh;
        
        % load base lfptrials
        files = dir(fullfile(inputfolder_lfpSKT, ['*_' basepd '_*.mat']));
        eval(['t_minmax_reach = t_minmax_reach_' basepd ';']);
        if segVFolder
            [lfptrials_base, fs_base, T_chnsarea_base] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
        else
            [lfptrials_base, fs_base, T_chnsarea_base] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
        end
        clear files
        
        % extract lfptrials_base of chns of AOI
        mask_chnsOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea_base.brainarea);
        lfptrials_base = lfptrials_base(mask_chnsOfI, :, :);
        T_chnsarea_base = T_chnsarea_base(mask_chnsOfI, :);
        clear mask_chnsOfI
        
        for ci = ci_str : ci_end
            pdcond = cond_cell{ci};
            if strcmpi(basepd, pdcond)
                continue;
            end
            
            [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond);
            disp([codefilename ': ' animal '-' pdcond '-' event])
            
            % saved file ciCohChangesfile
            ciCohChangesfile = fullfile(savematfolder, [ciCohChangesfile_prefix '_' basepd '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz' '_' pdcond '_' event '_align2' align2name '.mat']);
            
            % generate ciCohPhasefile if not exist
            if (~exist(ciCohChangesfile, 'file'))
                
                % load pdcond ciCoh
                load(fullfile(inputfolder_ciCohSKT, [animal ' ciCohPhasefile' '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz_' pdcond '_' event '_align2' align2name '.mat']), 'ciCoh', 'T_chnsarea', 'f_selected', 'ntrials');
                ciCoh_pd = ciCoh;
                clear ciCoh
                
                % calc ciCOHChanges
                ciCohChanges = ciCoh_pd - ciCoh_base;
                
                % load lfptrials of pdcond
                eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
                files = dir(fullfile(inputfolder_lfpSKT, ['*_' pdcond '_*.mat']));
                if segVFolder
                    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
                else
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
                end
                
                
                
                % extract lfptrials of chnOfI
                mask_chnsOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                lfptrials = lfptrials(mask_chnsOfI, :, :);
                T_chnsarea = T_chnsarea(mask_chnsOfI, :);
                clear mask_chnsOfI
                
                if ~isequal(fs, fs_base) || ~isequal(T_chnsarea, T_chnsarea_base)
                    disp([event '-' pdcond 'fs or T_chnsarea not equal to base' ])
                    continue;
                end
                
                % save before psedo Test
                save(ciCohChangesfile, 'ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected');
                clear('ciCohChanges', 'T_chnsarea', 'f_selected', 'ntrials');
                
                % psedo Test
                psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile);
                
                clear files lfptrials fs 
            end
            
            load(ciCohChangesfile, 'ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected', 'psedoiCohChanges');
            
            if ~exist('psedoiCohChanges','var') || size(psedoiCohChanges, 4) < shuffleN_psedoTest
                eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
                files = dir(fullfile(inputfolder_lfpSKT, ['*_' pdcond '_*.mat']));
                if segVFolder
                    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach);
                else
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach); % lfptrials: nchns * ntemp * ntrials
                end
                % extract data of chns of AOI
                mask_chnsOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                lfptrials = lfptrials(mask_chnsOfI, :, :);
                T_chnsarea = T_chnsarea(mask_chnsOfI, :);
                clear mask_chnsOfI
                
                %  extract and save psedociCohs
                psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile);
                
                clear lfptrials fs_lfpSKT T_chnsarea_lfpSKT files t_minmax_reach
                
                load(ciCohChangesfile, 'psedoiCohChanges');
            end
            nshuffle = size(psedoiCohChanges, 4);
            
            % extract sigciCoh
            [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
            
            % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
            [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
            [ciCohChanges_flatten_used, chnPairNames_used]= ciCoh_Used(chnPairNames, sigciCohChanges_flatten, []);
            
            
            % plot and save ciCohChanges Histogram image
            titlename = [animal ' FC Changes '  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            plot_ciCohHistogram(ciCohChanges_flatten_used, chnPairNames_used, round(f_selected), titlename, ...
                                'histClim', [-0.5 0.5], 'fig_width', 1000, 'fig_height', 250);
            saveimgname = [ciCohChangesfile_prefix '_' basepd  '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz' '_' event '_' pdcond '_align2' char(align2) '.' image_type];
            saveas(gcf, fullfile(saveimgfolder, saveimgname), image_type);
            clear titlename  saveimgname
            close all
            
            clear align2 t_AOI align2name ciCohChangesfile
            clear('ciCohChanges', 'T_chnsarea', 'ntrials', 'f_selected', 'psedoiCohChanges');
        end
        
    end
    clear f_AOI
end

