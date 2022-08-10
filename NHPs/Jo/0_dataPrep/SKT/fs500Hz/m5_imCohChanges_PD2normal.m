function m5_imCohChanges_PD2normal(varargin)
% combine mild and moderate state
%
%   Example Usage:
%           m5_imCohChanges_PD2normal('newRun', true)
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
inputfolder_ciCoh_normal = fullfile(codecorresParentfolder, 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');
inputfolder_lfp_normal = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');
inputfolder_ciCoh_PD = fullfile(codecorresParentfolder, 'm4_imCoh_UsingFFT_combinePD');



%% save setup
savefolder = codecorresfolder;

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

ciCohChangesfile_prefix =[animal '-ciCohChanges_PD2Normal_'];



%% Code start here
EventPhases = SKT_eventPhases_extract(animal, 'codesavefolder', savecodefolder);
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
for ei = 1 : length(EventPhases)
    event = EventPhases{ei};

    disp([codefilename ' ' animal '-' event])

    [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, 'codesavefolder', savecodefolder);
    ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix  '_PD2Normal_' event '_align2' align2name '.mat']);

    if(~exist(ciCohChangesfile, 'file') || newRun)
        % load PD related variables
        ciCohfile_PD = fullfile(inputfolder_ciCoh_PD, [animal '_ciCoh__PD_'  event '_align2' align2name '.mat']);
        load(ciCohfile_PD, 'ciCohs', 'fs', 'f_selected', 'lfptrials', 'psedociCohs', 'T_chnsarea')
        ciCohs_PD = ciCohs;
        fs_PD = fs;
        f_selected_PD = f_selected;
        lfptrials_PD = lfptrials;
        psedociCohs_PD = psedociCohs;
        T_chnsarea_PD = T_chnsarea;
        clear('ciCohs', 'fs', 'f_selected', 'lfptrials', 'psedociCohs', 'T_chnsarea')

        
        % load normal related variables
        ciCohfile_normal = fullfile(inputfolder_ciCoh_normal, [animal ' ciCohPhasefile_normal_' event '_align2' align2name '.mat']);
        load(ciCohfile_normal, 'ciCoh','T_chnsarea', 'psedociCohs', 'f_selected');
        chnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
        ciCohs_normal = ciCoh(chnsOfI_mask, chnsOfI_mask, :);
        psedociCohs_normal = psedociCohs(chnsOfI_mask, chnsOfI_mask, :, :);
        T_chnsarea_normal = T_chnsarea(chnsOfI_mask, :);
        f_selected_normal = f_selected;
        files_normal = dir(fullfile(inputfolder_lfp_normal, ['*_normal_*.mat']));
        if(isempty(files_normal))
            disp("normal lfp files are empty")
        end
        [lfptrials_normal, ~, ~] = lfp_goodTrials_align2(files_normal, align2, t_AOI, [0.5 1]);
        lfptrials_normal = lfptrials_normal(chnsOfI_mask, :, :);
        clear chnsOfI_mask files_normal
        clear('ciCoh','T_chnsarea', 'psedociCohs', 'f_selected')


        % save exist variables
        ciCohs.PD = ciCohs_PD;
        ciCohs.normal = ciCohs_normal;
        lfptrials.PD = lfptrials_PD;
        lfptrials.normal = lfptrials_normal;
        psedociCohs.PD = psedociCohs_PD;
        psedociCohs.normal = psedociCohs_normal;
        if(isequal(f_selected_PD, f_selected_normal) && isequal(T_chnsarea_PD, T_chnsarea_normal))
            fs = fs_PD;
            f_selected = f_selected_PD;
            T_chnsarea = T_chnsarea_PD;
        end
        
        save(ciCohChangesfile, 'ciCohs', 'psedociCohs', 'lfptrials', 'fs', 'f_selected', 'T_chnsarea') 

        % save ciCohChanges
        ciCohChanges = ciCohs.PD - ciCohs.normal;
        save(ciCohChangesfile, 'ciCohChanges', '-append')
        
        clear ciCohs_PD  lfptrials_PD  psedociCohs_PD fs_PD f_selected_PD   T_chnsarea_PD 
        clear ciCohs_normal lfptrials_normal psedociCohs_normal f_selected_normal T_chnsarea_normal
        clear('ciCohs', 'psedociCohs', 'lfptrials', 'fs', 'f_selected', 'T_chnsarea')
    end

    % psedociCohChanges
    load(ciCohChangesfile, 'lfptrials', 'fs') 
    psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials.PD, lfptrials.normal, fs, f_AOI, ciCohChangesfile);
    clear('lfptrials', 'fs');

    clear event align2 t_AOI align2name ciCohChangesfile
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