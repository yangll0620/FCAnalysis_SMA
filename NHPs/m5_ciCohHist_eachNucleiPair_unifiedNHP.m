function m5_ciCohHist_eachNucleiPair_unifiedNHP(animal, varargin)
%
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

cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();

% parse params
p = inputParser;
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', length(EventPhases), @isscalar);
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', length(cond_cell), @isscalar);
parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;

% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
segVFolder = false;
if strcmpi(animal, 'Kitty')
    segVFolder = true;
end
    
if segVFolder  
    SKTSubfolder = 'SKT_SegV';
    SKTDataSubfolder = 'm2_segSKTData_SelectTrials_goodReach';
else
    SKTSubfolder = 'SKT';
    SKTDataSubfolder = 'm2_SKTData_SelectTrials';
end
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, codefilename);
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

%% Input setup
inputfolder = fullfile(codecorresParentfolder, 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');

image_type = 'tif';

f_AOI = [8 40];

histClim = [0 1];

runCicohHist = true;


%% Code start here

unwanted_DBS = unwanted_DBS_extract(animal,'codesavefolder', savecodefolder);
noisy_chns = noisy_chns_extract(animal,'codesavefolder', savecodefolder);
notAOI_chns = notInterested_chns_extract(animal, 'codesavefolder', savecodefolder);
removed_chns = [unwanted_DBS noisy_chns notAOI_chns];
clear unwanted_DBS noisy_chns

for ei = ei_str: ei_end
    event = EventPhases{ei};
    
    subphsavefolder = fullfile(phsubfolder, event);
    if ~exist(subphsavefolder, 'dir')
        mkdir(subphsavefolder);
    end
    
    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};
        
        disp([codefilename ' ' animal '-' event '-' pdcond])
        
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, 'codesavefolder', savecodefolder);
        
        % load(and extract) ciCohPhasefile
        ciCohPhasefile = fullfile(inputfolder, [ciCohPhasefile_prefix  '_' pdcond '_' event '_align2' align2name '.mat']);

        load(ciCohPhasefile, 'ciCoh', 'deltaphis_allChnsTrials','T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');

        % extract sigciCoh
        [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh, 'codesavefolder', savecodefolder);
        

        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [sigciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCoh, deltaphis_allChnsTrials, T_chnsarea, 'codesavefolder', savecodefolder);
        [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, sigciCoh_flatten, deltaphis_flatten, removed_chns, 'codesavefolder', savecodefolder);
        
        
        % plot and save ciCoh Histogram image
        if runCicohHist
            nshuffle = size(psedociCohs,4);
            
            % M1-STN pair
            pairClus = 'M1STN';
            mask = cellfun(@(x) contains(x, 'M1')&&contains(x, 'stn'), chnPairNames_used);
            titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s ' pairClus ', align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            plot_ciCohHistogram(ciCoh_flatten_used(mask, :), chnPairNames_used(mask), f_selected, titlename, histClim, 'codesavefolder', savecodefolder);
            saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '_' pairClus '.' image_type];
            saveas(gcf, fullfile(savefolder, saveimgname), image_type);
            clear pairClus mask titlename  saveimgname
            close(gcf)
            
            % M1-GP pair
            pairClus = 'M1GP';
            mask = cellfun(@(x) contains(x, 'M1')&&contains(x, 'gp'), chnPairNames_used);
            titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s ' pairClus ', align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            plot_ciCohHistogram(ciCoh_flatten_used(mask, :), chnPairNames_used(mask), f_selected, titlename, histClim, 'codesavefolder', savecodefolder);
            saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '_' pairClus '.' image_type];
            saveas(gcf, fullfile(savefolder, saveimgname), image_type);
            clear pairClus mask titlename  saveimgname
            close(gcf)
            
            
            % STN-GP pair
            pairClus = 'STNGP';
            mask = cellfun(@(x) contains(x, 'stn')&&contains(x, 'gp'), chnPairNames_used);
            titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s ' pairClus ', align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            plot_ciCohHistogram(ciCoh_flatten_used(mask, :), chnPairNames_used(mask), f_selected, titlename, histClim, 'codesavefolder', savecodefolder);
            saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '_' pairClus '.' image_type];
            saveas(gcf, fullfile(savefolder, saveimgname), image_type);
            clear pairClus mask titlename  saveimgname
            close(gcf)
        end
        
 
         
        clear pdcond subpdsavefolder align2 t_AOI align2name ciCohPhasefile
        clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
    end
end
