function m5_imCoh_Freeze(varargin)
% calculate and plot cicoh of freeze, here freeze contains init, reach and
% mani freeze together
%
%   Example Usage:
%           m5_imCohChanges_reachFreeze2ReachPhases('shuffleN_psedoTest', 500, 'newRun', true, 'plotCiCoh', true)
%   
%   Input:
%       Name-Value:
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500
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
addParameter(p, 'plotCiCoh', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
plotCiCoh = p.Results.plotCiCoh;



% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%%  input setup

lfpfile = fullfile(codecorresParentfolder, 'm4_fs500Hz_FreezeSegs_extract', 'Kitty-FreezeSegs-ReachTrials-FreezeSegs.mat');

f_AOI = [8 40];


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code start here
ciCoh_file = fullfile(savefolder, 'ciCohs-Freeze.mat');


if (~exist(ciCoh_file, 'file') || newRun)
    load(lfpfile, 'lfpsegs_Freeze', 'fs', 'T_chnsarea');

    % extract all lfpsegs from all freeze types
    lfpsegs_InitFreeze = cat(3, lfpsegs_Freeze.InitFreeze.earlyFreeze, lfpsegs_Freeze.InitFreeze.middleFreeze, lfpsegs_Freeze.InitFreeze.lateFreeze); 
    lfpsegs_ReachFreeze = cat(3, lfpsegs_Freeze.ReachFreeze.earlyFreeze, lfpsegs_Freeze.ReachFreeze.middleFreeze, lfpsegs_Freeze.ReachFreeze.lateFreeze); 
    lfpsegs_ManiFreeze = cat(3, lfpsegs_Freeze.ManiFreeze.earlyFreeze, lfpsegs_Freeze.ManiFreeze.middleFreeze, lfpsegs_Freeze.ManiFreeze.lateFreeze); 
    lfpsegs = cat(3, lfpsegs_InitFreeze, lfpsegs_ReachFreeze, lfpsegs_ManiFreeze);
    clear lfpsegs_InitFreeze lfpsegs_ReachFreeze lfpsegs_ManiFreeze

    % select 1000
    inds = randsample(size(lfpsegs, 3),100);
    lfpsegs  = lfpsegs(:, :, inds);
    clear inds

    [~, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfpsegs, fs, f_AOI);

    ciCohs.freeze = ciCoh;
    clear ciCoh

    save(ciCoh_file, 'lfpsegs',  'fs', 'f_AOI','T_chnsarea');
    save(ciCoh_file, 'ciCohs', 'f_selected',  '-append');

    clear('lfpsegs_Freeze', 'fs', 'T_chnsarea');
    clear('ciCohs', 'f_selected');
end

% psedoCicoh and psedoCicohChanges
load(ciCoh_file, 'lfpsegs', 'fs', 'f_AOI')
psedoFreezeCiCoh_extract_save(shuffleN_psedoTest, lfpsegs, fs, f_AOI, ciCoh_file);
clear('lfpsegs', 'fs', 'f_AOI')



function psedoFreezeCiCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohfile)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohfile:
%             
%   
% save psedociCohs.(slowFastType).(ePhase): nchns * nchns * nf * nshuffle, saved to ciCohfile


nchns = size(lfptrials, 1);

load(ciCohfile, 'psedociCohs');
if ~exist('psedociCohs', 'var')
    psedociCohs = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs, 4) + 1;
end

disp('Freeze psedo test ')
for si = shuffi_str : suffi_end
    for chni = 1 : nchns -1
        lfp1 = squeeze(lfptrials(chni, :, :));
        for chnj = chni + 1 : nchns
            lfp2 = squeeze(lfptrials(chnj, :, :));
            [psedociCoh, ~] = psedo_ciCohSKT_FFT_NoAmp(lfp1, lfp2, fs, f_AOI);
            
            if chni == 1 && chnj == chni + 1
                nf = size(psedociCoh, 1);
                psedociCohs1Time = zeros(nchns, nchns, nf, 1);
                clear nf
            end
            
            psedociCohs1Time(chni, chnj, :, :) = psedociCoh;
            clear lfp2 psedociCoh
        end
        clear lfp1
    end
    psedociCohs = cat(4, psedociCohs, psedociCohs1Time);
    clear psedociCohs1Time

    if(mod(si,100) == 0)
        disp(['Freeze psedo test ' num2str(si) ' times'])
        save(ciCohfile, 'psedociCohs', '-append');
    end
end

