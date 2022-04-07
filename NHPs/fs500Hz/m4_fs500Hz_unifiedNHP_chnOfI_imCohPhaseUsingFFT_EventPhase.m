function m4_fs500Hz_SKTData_imCohPhaseUsingFFT_EventPhase(animal, varargin)
%
%   Input:
%       animal
%
%       Name-Value: 
%           ei_str - event start index
%           ei_end - event end index
%           ci_str - condition start index
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
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));


cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();

% parse params
p = inputParser;
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', length(EventPhases), @isscalar);
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', length(cond_cell), @isscalar);
addParameter(p, 'ishighfreq', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'runCicohHist', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'runRosePlot', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));



parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
ishighfreq = p.Results.ishighfreq;
runCicohHist = p.Results.runCicohHist;
runRosePlot = p.Results.runRosePlot;
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
    SKTDataSubfolder = 'm2_segSKTData_SelectTrials_goodReach_fs1000Hz';
else
    SKTSubfolder = 'SKT';
    SKTDataSubfolder = 'm2_SKTData_SelectTrials';
end
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, 'fs500Hz', codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, 'dir')
    rmdir(savecodefolder,'s');
end
copyfile2folder(codefilepath, savecodefolder);

phsubfolder = fullfile(savefolder, 'phases');
if ~exist(phsubfolder, 'dir')
    mkdir(phsubfolder);
end


ciCohPhasefile_prefix =[animal ' ciCohPhasefile'];

%%  input setup

inputfolder = fullfile(codecorresParentfolder, SKTDataSubfolder);
image_type = 'tif';

if runCicohHist
    histClim = [0 1];
end
if runRosePlot
    roseRLim = [0 0.3];
end


%% Code start here
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal, 'codesavefolder', savecodefolder);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', savecodefolder);

f_AOIs = freqsOfInterest_extract(animal, 'ishighfreq', ishighfreq, 'codesavefolder', fullfile(savefolder, 'code'));
disp(['f_AOIs = ' num2str(f_AOIs(1)) ' ' num2str(f_AOIs(2)) 'Hz'] )
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
for fi = 1: size(f_AOIs, 1)
    f_AOI = f_AOIs(fi, :);
    
    for ei = ei_str: ei_end
        event = EventPhases{ei};

        subphsavefolder = fullfile(phsubfolder, event);
        if ~exist(subphsavefolder, 'dir')
            mkdir(subphsavefolder);
        end

        for ci = ci_str : ci_end
            pdcond = cond_cell{ci};

            disp([codefilename ': ' animal ' ' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz ' event '-' pdcond])

            [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, 'codesavefolder', savecodefolder);

            % load(and extract) ciCohPhasefile
            ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz_' pdcond '_' event '_align2' align2name '.mat']);
            if(~exist(ciCohPhasefile, 'file'))

                files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
                if isempty(files)
                    continue;
                end

                eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
                if segVFolder
                    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder',savecodefolder);
                else
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder); % lfptrials: nchns * ntemp * ntrials
                end

                % extract data of chns of AOI
                mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                lfptrials = lfptrials(mask_chnOfI, :, :);
                T_chnsarea = T_chnsarea(mask_chnOfI, :);

                %  extract and save deltaphis_allChnsTrials and cicoh
                [deltaphis_allChnsTrials, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfptrials, fs, f_AOI, 'codesavefolder', savecodefolder);


                ntrials = size(lfptrials, 3);
                save(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');

                %  extract and save psedociCohs
                psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);


                clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
                clear t_minmax_reach files lfptrials
            end
            load(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
            if ~exist('psedociCohs','var')
                eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
                files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
                if segVFolder
                    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder);
                else
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder); % lfptrials: nchns * ntemp * ntrials
                end

                % extract data of chns of AOI
                mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                lfptrials = lfptrials(mask_chnOfI, :, :);
                T_chnsarea = T_chnsarea(mask_chnOfI, :);

                psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
                load(ciCohPhasefile, 'psedociCohs');
                clear t_minmax_reach lfptrials
            end

            nshuffle = size(psedociCohs, 4);
            if(nshuffle < shuffleN_psedoTest)
                eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
                files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
                if segVFolder
                    [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder);
                else
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, [t_AOI(1) t_AOI(2)], t_minmax_reach, 'codesavefolder', savecodefolder); % lfptrials: nchns * ntemp * ntrials
                end

                % extract data of chns of AOI
                mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                lfptrials = lfptrials(mask_chnOfI, :, :);
                T_chnsarea = T_chnsarea(mask_chnOfI, :);


                psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
                load(ciCohPhasefile, 'psedociCohs');
                clear t_minmax_reach lfptrials
            end
            nshuffle = size(psedociCohs, 4);

            % extract sigciCoh
            [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh, 'codesavefolder', savecodefolder);


            % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
            [sigciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCoh, deltaphis_allChnsTrials, T_chnsarea, 'codesavefolder', savecodefolder);
            [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, sigciCoh_flatten, deltaphis_flatten, [], 'codesavefolder', savecodefolder);


            % plot and save ciCoh Histogram image
            if runCicohHist
                titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
                plot_ciCohHistogram(ciCoh_flatten_used, chnPairNames_used, f_selected, titlename, 'histClim', histClim, 'codesavefolder', savecodefolder);
                saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type];
                saveas(gcf, fullfile(savefolder, saveimgname), image_type);
                clear titlename  saveimgname
            end


            % rose histogram of deltaphis_allChnsTrials
            if runRosePlot
                titlename_prefix = [animal '-'  pdcond '-'  event];
                subtitlename = [event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s, align2 = ' char(align2) ', ntrials = ' num2str(ntrials)];
                savefile_prefix = [animal 'trialPhaseDiff'];
                savefile_suffix = [event '_' pdcond '_align2' char(align2)];
                plotsave_deltaphirose(deltaphis_flatten_used, ciCoh_flatten_used, chnPairNames_used, f_selected, titlename_prefix, subtitlename, subphsavefolder, savefile_prefix, savefile_suffix, image_type,...
                    'codesavefolder', savecodefolder, 'roseRLim', roseRLim);
                clear titlename_prefix subtitlename savefile_prefix savefile_suffix
            end


            clear pdcond subpdsavefolder align2 t_AOI align2name ciCohPhasefile
            clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
            
            close all
        end
    end
    
    clear f_AOI
end
