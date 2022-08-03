function m6_imCohChangesChanges_reachFreeze2Reach(varargin)
% calculate ciCohChangesofChanges for example: 
%   ciCohChanges.b_earlyReach.ReachFreeze.earlyFreeze = (ciCohs.ReachFreeze.earlyFreeze - ciCoh_earlyReach);
%
%
%   ciCohChangesChanges_bearlyReach_middle2earlyFreeze = 
%   ciCohChanges.b_earlyReach.ReachFreeze.middleFreeze - ciCohChanges.b_earlyReach.ReachFreeze.earlyFreeze;
%
%   Example Usage:
%           m6_imCohChangesChanges_reachFreeze2Reach();
%           m6_imCohChangesChanges_reachFreeze2Reach('shuffleN_psedoTest', 500, 'newRun', true)
%           m6_imCohChangesChanges_reachFreeze2Reach('plotCiCohChanges', false, 'plotCiCoh', true)
%   
%   Input:
%       Name-Value:
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500
%           plotCiCohChanges - plot ciCohChanges true (default) or false
%           plotCiCoh - plot ciCoh true (default) or false

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
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plotCiCohChanges', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'plotCiCoh', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
plotCiCohChanges = p.Results.plotCiCohChanges;
plotCiCoh = p.Results.plotCiCoh;



% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


% animal
animal = animal_extract(codecorresfolder);


%%  input setup
ciCoh_ChangesFile = fullfile(codecorresParentfolder, 'm5_imCohChanges_reachFreeze2ReachPhases', 'ciCohs-Changes-reachFreeze2reachPhases.mat');



%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code start here
ciCoh_ChangesChanges_file = fullfile(savefolder, 'ciCohs-ChangesChanges-reachFreezeAlongTime.mat');
if (~exist(ciCoh_ChangesChanges_file, 'file') || newRun)
    load(ciCoh_ChangesFile, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'T_chnsarea', 'f_selected', ...
        'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'f_AOI');


    save(ciCoh_ChangesChanges_file, 'ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'T_chnsarea', 'f_selected', ...
        'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'f_AOI');


    % calculate ciCohChangesChanges
    eBasePhases = fieldnames(ciCohChanges);
    for ei = 1 : length(eBasePhases)
        eBasePhase = eBasePhases{ei};

        reachfreezeTypes = fieldnames(ciCohChanges.(eBasePhase).ReachFreeze);
        for fri = 1 : length(reachfreezeTypes)-1
            subfreezeType_base = reachfreezeTypes{fri};

            ciCohChanges_base = ciCohChanges.(eBasePhase).ReachFreeze.(subfreezeType_base);

            for frj = fri+1 : length(reachfreezeTypes)
                subfreezeType_comp = reachfreezeTypes{frj};
                ciCohChanges_comp = ciCohChanges.(eBasePhase).ReachFreeze.(subfreezeType_comp);

                ciCohChangesChanges.(eBasePhase).ReachFreeze.([subfreezeType_comp '2' subfreezeType_base]) = ciCohChanges_comp - ciCohChanges_base;

                clear subfreezeType_comp ciCohChanges_comp
            end

            clear subfreezeType_base ciCohChanges_base
        end

        clear eBasePhase reachfreezeTypes
    end

    save(ciCoh_ChangesChanges_file, 'ciCohChangesChanges', '-append');
    clear ciCohChangesChanges

    clear('ciCohs', 'psedociCohs', 'ciCohChanges', 'psedociCohChanges', 'T_chnsarea', 'f_selected', ...
        'lfpsegs_Freeze', 'lfptrials_Reach', 'seg_tseg', 'fs', 'f_AOI');
end


% psedo ciCohChangesChanges
load(ciCoh_ChangesChanges_file, 'ciCohChangesChanges', 'lfptrials_Reach', 'lfpsegs_Freeze', 'fs', 'f_AOI');
eBasePhases = fieldnames(ciCohChangesChanges);
reachfreezeTypes = fieldnames(ciCohChangesChanges.(eBasePhases{1}).ReachFreeze);
for ebi = 1 : length(eBasePhases)
    eBasePhase = eBasePhases{ebi};
    lfptrials_baseReach = lfptrials_Reach.(eBasePhase);
    for fri = 1 : length(reachfreezeTypes)-1
        subfreezeType_base = reachfreezeTypes{fri};
        lfpsegs_freeze_base = lfpsegs_Freeze.ReachFreeze.(subfreezeType_base);
        for frj = fri+1 : length(reachfreezeTypes)
            subfreezeType_comp = reachfreezeTypes{frj};
            lfpsegs_freeze_comp = lfpsegs_Freeze.ReachFreeze.(subfreezeType_comp);

            psedoCiCohChangesChanges_extract_save(shuffleN_psedoTest, lfpsegs_freeze_comp, lfpsegs_freeze_base, lfptrials_baseReach, fs, f_AOI, ciCoh_ChangesChanges_file, ...
                eBasePhase, 'ReachFreeze', [subfreezeType_comp '2' subfreezeType_base]);

            clear subfreezeType_comp lfpsegs_freeze_comp
        end

        clear subfreezeType_base lfpsegs_freeze_base
    end

    clear eBasePhase lfptrials_baseReach
end
clear('ciCohChangesChanges', 'lfptrials_Reach', 'lfpsegs_Freeze', 'fs', 'f_AOI');




function psedoCiCohChangesChanges_extract_save(suffi_end, lfptrials, lfptrials_base, lfptrials_beMinuend,fs, f_AOI, ciCohChangesChangesfile,  eBasePhase, freezeType, subfreezeType2base)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohPhasefile:
%       
%   
% save psedoiCohChanges: nchns * nchns * nf * nshuffle, saved to ciCohChangesfile



load(ciCohChangesChangesfile, 'psedociCohChangesChanges');

if(~exist('psedociCohChanges', 'var'))
    psedociCohChangesChanges = struct();
end

if ~isfield(psedociCohChangesChanges, eBasePhase) || ~isfield(psedociCohChangesChanges.(eBasePhase), freezeType) || ~isfield(psedociCohChangesChanges.(eBasePhase).(freezeType), subfreezeType2base)
    psedociCohChangesChanges.(eBasePhase).(freezeType).(subfreezeType2base) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChangesChanges.(eBasePhase).(freezeType).(subfreezeType2base), 4) + 1;
end

[~, ciCoh_beMinuend, ~] = ciCoh_trialDeltaPhi(lfptrials_beMinuend, fs, f_AOI);

lfp_combined = cat(3, lfptrials_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_Base = lfp_combined(:, :, masksBase);
    psedolfp_comp = lfp_combined(:, :, ~masksBase);
    
    [~, psedoiCoh_comp, ~] = ciCoh_trialDeltaPhi(psedolfp_comp, fs, f_AOI);
    [~, psedoiCoh_Base, ~] = ciCoh_trialDeltaPhi(psedolfp_Base, fs, f_AOI);

    psediChanges = (psedolfp_comp - ciCoh_beMinuend) - (psedolfp_Base - ciCoh_beMinuend);
    
    psedociCohChangesChanges.(eBasePhase).(freezeType).(subfreezeType2base) = cat(4, psedociCohChangesChanges.(eBasePhase).(freezeType).(subfreezeType2base), psedoiCoh_comp - psedoiCoh_Base);
    
    if(mod(si, 100) == 0)
        disp([eBasePhase '-' freezeType '-' subfreezeType2base ' pesdo ciCoh Changes test ' num2str(si)])
        save(ciCohChangesChangesfile, 'psedociCohChangesChanges', '-append');
    end
    
    clear randomSKTInds randomRestInds psedolfp_comp psedoiCoh_rest
    clear psedoiCoh_comp psedoiCoh_rest
end

