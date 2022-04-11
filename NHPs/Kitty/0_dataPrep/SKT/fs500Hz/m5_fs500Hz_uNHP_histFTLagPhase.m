function m5_fs500Hz_uNHP_histFTLagPhase(animal, varargin)
% plot cicoh Histogram, frequency time lag and phase of interested channels
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
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', [], @isscalar);
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', [], @isscalar);
addParameter(p, 'runCicohHist', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'runRosePlot', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end

parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
runCicohHist = p.Results.runCicohHist;
runRosePlot = p.Results.runRosePlot;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
disp('p.Results =  ' )
p.Results


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

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
inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');

image_type = 'tif';
f_AOI = [8 40];

if runCicohHist
    histClim = [0 1];
    
    histsavefolder = fullfile(savefolder, 'ciCohHist');
    if ~exist(histsavefolder, 'dir')
        mkdir(histsavefolder)
    end
end
if runRosePlot
    roseRLim = [0 0.3];
    
    savefile_prefix = [animal 'trialPhaseDiff'];
end


%% Code start here
cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();
if isempty(ei_end)
    ei_end = length(EventPhases);
end
if isempty(ci_end)
    ci_end = length(cond_cell);
end
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

for ei = ei_str: ei_end
    event = EventPhases{ei};
    
    subphsavefolder = fullfile(phsubfolder, event);
    if ~exist(subphsavefolder, 'dir')
        mkdir(subphsavefolder);
    end
    
    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};

        disp([ animal ' ' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz ' event '-' pdcond])

        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, 'codesavefolder', savecodefolder);

        % load(and extract) ciCohPhasefile
        ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz_' pdcond '_' event '_align2' align2name '.mat']);
        
        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ciCohPhasefile, 'file') || newRun)

            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end

            % lfpseg align2
            [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2,t_AOI, 'codesavefolder', savecodefolder);

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


            clear files lfptrials fs T_chnsarea mask_chnOfI
            clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
        end
        
        %%% ----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
        pseciCohsVar = whos('-file',ciCohPhasefile, 'psedociCohs');
        if isempty(pseciCohsVar) || pseciCohsVar.size(4)< shuffleN_psedoTest
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end
            
            % lfpseg align2
            [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2,t_AOI, 'codesavefolder', savecodefolder);

            % extract data of chns of AOI
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            lfptrials = lfptrials(mask_chnOfI, :, :);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);

            %  extract and save psedociCohs
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
            
            clear files lfptrials fs T_chnsarea mask_chnOfI
        end
        clear pseciCohsVar


        %%% -- plot section --- %%%
        load(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
        
        % extract sigciCoh
        [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh, 'codesavefolder', savecodefolder);


        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [sigciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCoh, deltaphis_allChnsTrials, T_chnsarea, 'codesavefolder', savecodefolder);
        [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, sigciCoh_flatten, deltaphis_flatten, [], 'codesavefolder', savecodefolder);


        % plot and save ciCoh Histogram image
        if runCicohHist
            nshuffle = size(psedociCohs, 4);
            titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            
            plot_ciCohHistogram(ciCoh_flatten_used, chnPairNames_used, f_selected, titlename, 'codesavefolder', savecodefolder);
            
            saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type];
            saveas(gcf, fullfile(histsavefolder, saveimgname), image_type);
            close gcf
            clear titlename  saveimgname nshuffle
        end


        % rose histogram of deltaphis_allChnsTrials
        if runRosePlot
            titlename_prefix = [animal '-'  pdcond '-'  event];
            subtitlename = [event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s, align2 = ' char(align2) ', ntrials = ' num2str(ntrials)];
            savefile_suffix = [event '_' pdcond '_align2' char(align2)];
            rosePlotsavefolder = fullfile(savefolder, 'rosePlot', ephase);
            if ~exist(rosePlotsavefolder, 'dir')
                mkdir(rosePlotsavefolder);
            end
            
            plotsave_deltaphirose(deltaphis_flatten_used, ciCoh_flatten_used, chnPairNames_used, f_selected, titlename_prefix, subtitlename, rosePlotsavefolder, savefile_prefix, savefile_suffix, image_type,...
                'codesavefolder', savecodefolder, 'roseRLim', roseRLim);
            close gcf
            clear titlename_prefix subtitlename savefile_prefix savefile_suffix
        end


        clear pdcond subpdsavefolder align2 t_AOI align2name ciCohPhasefile
        clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');

        close all
    end
end