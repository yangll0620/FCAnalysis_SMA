function m4_imCoh_UsingFFT_combinePD(varargin)
%
%   Example Usage:
%           m4_imCohPhaseUsingFFT_EventPhase_unifiedNHP('Jo', 'ei_str', 1, 'ci_str', 1)
%

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

% parse parameter
p = inputParser;
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;

% pipelinefolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%  animal 
animal = animal_extract(codecorresfolder);


%%  input setup
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, 'dir')
    rmdir(savecodefolder,'s');
end
copyfile2folder(codefilepath, savecodefolder);


ciCohfile_prefix =[animal '_ciCoh_'];

f_AOI = [8 40];



%% Code start here
t_minmax_reach_PD = [0.6 1];


chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
EventPhases = SKT_eventPhases_extract(animal);
pdConds = {'mild', 'moderate'};
for ei = 1: length(EventPhases)
    event = EventPhases{ei};

    disp([codefilename ' ' animal '-' event])

    [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, 'codesavefolder', savecodefolder);
    ciCohfile = fullfile(savefolder, [ciCohfile_prefix  '_PD_' event '_align2' align2name '.mat']);


    if(~exist(ciCohfile, 'file') || newRun)
       
        % extract combined pd lfp
        pdfiles = [];
        for ci = 1 : length(pdConds)
            pdcond = pdConds{ci};
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            pdfiles = [pdfiles; files];
            clear pdcond files
        end
        [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(pdfiles, align2, t_AOI, t_minmax_reach_PD);
        clear pdfiles

        % remove unused chns
        chnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
        lfptrials = lfptrials(chnsOfI_mask, :, :);
        T_chnsarea = T_chnsarea(chnsOfI_mask, :);
        clear chnsOfI_mask

        % save
        save(ciCohfile, 'lfptrials', 'T_chnsarea', 'fs');
        [ciCohs, f_selected]= ciCoh_trial(lfptrials, fs, f_AOI);
        save(ciCohfile, 'ciCohs', 'f_selected', '-append');
        
        clear lfptrials T_chnsarea fs ciCohs f_selected
    end
    
    % psedo cicoh
    load(ciCohfile, 'lfptrials', 'fs')
    psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohfile);
    clear lfptrials fs

end


function psedociCoh_extract_save(suffi_end, lfptrials, fs, f_AOI, ciCohPhasefile)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohPhasefile:
%              
%   
% save psedociCohs: nchns * nchns * nf * nshuffle, saved to ciCohPhasefile



nchns = size(lfptrials, 1);

load(ciCohPhasefile, 'psedociCohs');
if(~exist('psedociCohs', 'var'))
    psedociCohs = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohs, 4) + 1;
end

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
        disp(['psedo test finished ' num2str(si) ' times'])
        save(ciCohPhasefile, 'psedociCohs', '-append');
    end
end


function [lfptrials, fs_lfp, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
%         Args:
%             align2: the event to be aligned (e.g Event.REACHONSET or uint 1-5)
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%             t_minmax_reach : min and max reach/return (s) for selecting trials (e.g [0.5 1])
%   
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


coli_align2 = uint32(align2);

coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    filename = files(filei).name;
    file = fullfile(files(filei).folder, filename);
    
    load(file, 'lfpdata', 'T_idxevent_lfp', 'goodTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');

    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    

    ntrials = size(lfpdata, 3);
    for tri = 1: ntrials
        
        % ignore trials marked with 0
        if ~goodTrials(tri)
            continue
        end
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_minmax_reach(1) 
            clear t_reach
            continue
        end


        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial(idx_reachonset_ma: idx_reach_ma, tri));
            idx_peakV_ma = idx + idx_reachonset_ma -1;
            t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
            
            if t_reachonset2peakV < abs(tdur_trial(1)) || t_peakV2reach < tdur_trial(2)
                clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
                clear t_reachonset2peakV  t_peakV2reach
                continue;
            end
            
            % extract trial with t_dur
            idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
            idx_time0 = idx_peakV_lfp;
            
            clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
            clear t_reachonset2peakV  t_peakV2reach
            clear idx_peakV_lfp
        else
            idx_time0 = T_idxevent_lfp{tri, coli_align2}; 
        end

               
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs_lfp) + idx_time0;
        idxdur(1) = idxdur(1) + 1;
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
       
        
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach t_return idxdur lfp_phase_1trial idx_time0
    end
    
    clear filename file ntrials
    clear('lfpdata', 'T_idxevent_lfp', 'goodTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
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