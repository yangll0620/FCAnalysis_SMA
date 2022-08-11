function m5_imCohChanges_PD_compEvent(varargin)
% combine mild and moderate state
%
%   Example Usage:
%           m5_imCohChanges_PD_compEvent('newRun', true)
%   
%   Input:
%       animal
%
%       Name-Value:
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500

codefilepath = mfilename('fullpath');
[~, codefilename]= fileparts(codefilepath);

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

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;


[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%  animal 
animal = animal_extract(codecorresfolder);


f_AOI = [8 40];

%%  input setup
inputfolder_ciCoh_PD = fullfile(codecorresParentfolder, 'm4_imCoh_UsingFFT_combinePD');



%% save setup
savefolder = codecorresfolder;

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

ciCohChangesfile_prefix =[animal '-ciCohChanges_PD_'];

%% Code start here

[~, tbl_compEvents]= compCondEvents_extract(animal, 'codesavefolder', savecodefolder);

for compEi = 1: height(tbl_compEvents)
    baseEvent = tbl_compEvents.baseEvent{compEi};
    compEvent = tbl_compEvents.compEvent{compEi};

    disp([codefilename ' ' animal '-b' baseEvent '-' compEvent])

    ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix  '_PD_b' baseEvent  '--' compEvent '.mat']);
    if(~exist(ciCohChangesfile, 'file') || newRun)
        

        % load base PD related variables
        [~, ~, align2name_base] = SKT_EventPhase_align2_tAOI_extract(baseEvent, animal, 'codesavefolder', savecodefolder);
        ciCohfile_base = fullfile(inputfolder_ciCoh_PD, [animal '_ciCoh__PD_'  baseEvent '_align2' align2name_base '.mat']);
        load(ciCohfile_base, 'ciCohs', 'fs', 'f_selected', 'lfptrials', 'psedociCohs', 'T_chnsarea')
        ciCohs_base = ciCohs;
        fs_base = fs;
        f_selected_base = f_selected;
        lfptrials_base = lfptrials;
        psedociCohs_base = psedociCohs;
        T_chnsarea_base = T_chnsarea;
        clear('ciCohs', 'fs', 'f_selected', 'lfptrials', 'psedociCohs', 'T_chnsarea')
        clear align2name_base


        % load base PD related variables
        [~, ~, align2name_comp] = SKT_EventPhase_align2_tAOI_extract(compEvent, animal, 'codesavefolder', savecodefolder);
        ciCohfile_comp = fullfile(inputfolder_ciCoh_PD, [animal '_ciCoh__PD_'  compEvent '_align2' align2name_comp '.mat']);
        load(ciCohfile_comp, 'ciCohs', 'fs', 'f_selected', 'lfptrials', 'psedociCohs', 'T_chnsarea')
        ciCohs_comp = ciCohs;
        fs_comp = fs;
        f_selected_comp = f_selected;
        lfptrials_comp = lfptrials;
        psedociCohs_comp = psedociCohs;
        T_chnsarea_comp = T_chnsarea;
        clear('ciCohs', 'fs', 'f_selected', 'lfptrials', 'psedociCohs', 'T_chnsarea')
        clear align2name_comp

        


        % save exist variables
        ciCohs.(baseEvent) = ciCohs_base;
        lfptrials.(baseEvent) = lfptrials_base;
        psedociCohs.(baseEvent) = psedociCohs_base;

        ciCohs.(compEvent) = ciCohs_comp;
        lfptrials.(compEvent) = lfptrials_comp;
        psedociCohs.(compEvent) = psedociCohs_comp;

        if(isequal(f_selected_base, f_selected_comp) && isequal(T_chnsarea_base, T_chnsarea_comp) && isequal(fs_base, fs_comp))
            fs = fs_base;
            f_selected = f_selected_base;
            T_chnsarea = T_chnsarea_base;
        end
        
        save(ciCohChangesfile, 'ciCohs', 'psedociCohs', 'lfptrials', 'fs', 'f_selected', 'T_chnsarea') 

        % save ciCohChanges
        ciCohChanges = ciCohs.(compEvent) - ciCohs.(baseEvent);
        save(ciCohChangesfile, 'ciCohChanges', '-append')
        
        clear ciCohs_base  lfptrials_base  psedociCohs_base fs_base f_selected_base   T_chnsarea_base 
        clear ciCohs_comp lfptrials_comp psedociCohs_comp f_selected_comp T_chnsarea_comp fs_comp
        clear('ciCohs', 'psedociCohs', 'lfptrials', 'fs', 'f_selected', 'T_chnsarea')
    end

    % psedociCohChanges
    load(ciCohChangesfile, 'lfptrials', 'fs') 
    psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials.(compEvent), lfptrials.(baseEvent), fs, f_AOI, ciCohChangesfile);
    clear('lfptrials', 'fs');

    clear event align2 t_AOI align2name ciCohChangesfile
end


function psedociCohChanges_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile)
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



load(ciCohChangesfile, 'psedociCohChanges');

if(~exist('psedoiCohChanges', 'var'))
    psedociCohChanges = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges, 4) + 1;
end

lfp_combined = cat(3, lfptrials_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_base = lfp_combined(:, :, masksBase);
    psedolfp_comp = lfp_combined(:, :, ~masksBase);
    

    [psedoiCoh_base, ~] = ciCoh_trial(psedolfp_base, fs, f_AOI);
    [psedoiCoh_comp, ~] = ciCoh_trial(psedolfp_comp, fs, f_AOI);
    
    psedociCohChanges = cat(4, psedociCohChanges, psedoiCoh_comp - psedoiCoh_base);
    
    if(mod(si, 100) == 0)
        disp(['psedo test finished ' num2str(si) ' times'])
        save(ciCohChangesfile, 'psedociCohChanges', '-append');
    end
    
    clear masksBase
    clear randomSKTInds psedolfp_comp psedolfp_base
    clear psedoiCoh_comp psedoiCoh_base
end



function [ciCoh, f_selected]= ciCoh_trial(lfptrials, fs, f_AOI)
%
%
%  Inputs:
%       lfptrials: nchns * ntemp * ntrials
%       fs:
%       f_AOI
%         
%   
%   Output:
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials
%       ciCoh:  nchns * nchns * nf
%       f_selected: nf * 1



% extract ciCoh
[ciCoh, f_selected] = ciCohSKTAllchns_FFT_NoAmp(lfptrials, fs, f_AOI);