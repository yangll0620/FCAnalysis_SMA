function m4_fs500Hz_uNHP_imCohChanges_compEvent(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_fs500Hz_uNHP_imCohChanges_compEvent('Kitty', 'compEi_str', 1, 'ci_str', 1)
%           m4_fs500Hz_uNHP_imCohChanges_compEvent('Kitty', 'compEi_str', 1, 'ci_str', 1, 'ci_end', 1, 'shuffleN_psedoTest', 500, 'newRun', true)
%   
%   Input:
%       animal
%
%       Name-Value:
%           compEi_str - compare events start index
%           compEi_end - compare events end index
%           ci_str - condition start index
%           ci_end - condition end index
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
addParameter(p, 'compEi_str', 1, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'compEi_end', [], @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', [], @isscalar);
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
compEi_str = p.Results.compEi_str;
compEi_end = p.Results.compEi_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
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
cond_cell = cond_cell_extract(animal);
if isempty(ci_end)
    ci_end = length(cond_cell);
end
[~, tbl_compEvents]= compCondEvents_extract(animal);
if isempty(compEi_end)
    compEi_end = height(tbl_compEvents);
end
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

for compEi = compEi_str: compEi_end
    baseEvent = tbl_compEvents.baseEvent{compEi};
    compEvent = tbl_compEvents.compEvent{compEi};
    
    if ~strcmpi(baseEvent, 'preMove')
        continue;
    end

    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};
        
        [align2_base, tdur_base, align2name_base] = SKT_EventPhase_align2_tAOI_extract(baseEvent, animal, pdcond, 'codesavefolder', savecodefolder);
        [align2_comp, tdur_comp, align2name_comp] = SKT_EventPhase_align2_tAOI_extract(compEvent, animal, pdcond, 'codesavefolder', savecodefolder);
        
        ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix '_' pdcond '_b' baseEvent '--' compEvent  '_align2' align2name_comp '.mat']);


        disp([animal '-' compEvent '-based ' baseEvent '-' pdcond])

        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ciCohChangesfile, 'file') || newRun)
            
            ciCohfile_base = fullfile(inputfolder_SKTciCoh, [inputfile_cicoh_prefix pdcond '_' baseEvent '_align2' align2name_base '.mat']);
            load(ciCohfile_base, 'ciCoh', 'psedociCohs', 'ntrials', 'f_selected', 'T_chnsarea');
            ciCoh_base = ciCoh; 
            psedociCohs_base = psedociCohs;
            clear ciCohfile_base ciCoh psedociCohs
            
            ciCohfile_comp = fullfile(inputfolder_SKTciCoh, [inputfile_cicoh_prefix pdcond '_' compEvent '_align2' align2name_comp '.mat']);
            load(ciCohfile_comp, 'ciCoh', 'psedociCohs');
            ciCoh_comp = ciCoh; 
            psedociCohs_comp = psedociCohs;
            clear ciCohfile_comp ciCoh psedociCohs
            
            
            if strcmpi(animal, 'Jo')
                mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);

                T_chnsarea = T_chnsarea(mask_chnOfI, :);
                ciCoh_base = ciCoh_base(mask_chnOfI, mask_chnOfI, :);
                ciCoh_comp = ciCoh_comp(mask_chnOfI, mask_chnOfI, :);
                psedociCohs_base = psedociCohs_base(mask_chnOfI, mask_chnOfI, :, :);
                psedociCohs_comp = psedociCohs_comp(mask_chnOfI, mask_chnOfI, :, :);

                clear mask_chnOfI
            end

            save(ciCohChangesfile, 'ciCoh_base', 'ciCoh_comp', 'psedociCohs_base', 'psedociCohs_comp', 'ntrials', 'f_selected', 'T_chnsarea')
            
            ciCohChanges = ciCoh_comp - ciCoh_base;
            save(ciCohChangesfile, 'ciCohChanges', '-append');

            clear ciCohChanges
            clear('ciCoh_base', 'ciCoh_comp', 'psedociCohs_base', 'psedociCohs_comp', 'ntrials', 'f_selected', 'T_chnsarea');
        end
            
        
        %%%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
        psedoVar = whos('-file',ciCohChangesfile, 'psedoiCohChanges');
        if isempty(psedoVar) || psedoVar.size(4)< shuffleN_psedoTest
            
            % extract lfptrials
            lfpfiles = dir(fullfile(inputfolder_SKTlfp, ['*_' pdcond '_*.mat']));
            if isempty(lfpfiles)
                clear lfpfiles
                continue;
            end
            
            if strcmpi(animal, 'Kitty')
                [lfptrials_base, fs, T_chnsarea] = lfpseg_selectedTrials_align2(lfpfiles, align2_base,tdur_base, 'codesavefolder', savecodefolder);
                [lfptrials_comp, ~, ~] = lfpseg_selectedTrials_align2(lfpfiles, align2_comp,tdur_comp, 'codesavefolder', savecodefolder);
            end
            if strcmpi(animal, 'Jo')
                [lfptrials_base, fs, T_chnsarea] = lfptrials_goodTrials_align2(lfpfiles, align2_base,tdur_base);
                [lfptrials_comp, ~, ~] = lfptrials_goodTrials_align2(lfpfiles, align2_comp,tdur_comp);
            end
            
            
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);
            
            
            lfptrials_base = lfptrials_base(mask_chnOfI, :, :);
            lfptrials_comp = lfptrials_comp(mask_chnOfI, :, :);
            
            psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials_comp, lfptrials_base, fs, f_AOI, ciCohChangesfile, 'codesavefolder', savecodefolder);
            
            clear lfpfiles lfptrials_base fs T_chnsarea mask_chnOfI lfptrials_base lfptrials_comp
        end
        clear psedoVar

        %%% -- plot and save section --- %%%
        load(ciCohChangesfile, 'ciCohChanges', 'psedoiCohChanges', 'f_selected', 'ntrials', 'T_chnsarea')
        
        % extract sigciCoh
        [sigciCohChanges]= sigciCoh_extract(psedoiCohChanges, ciCohChanges);
        [sigciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
        
        % plot and save ciCoh Histogram image
        titlename = [imgtitle_prefix ':' pdcond '-' compEvent '/' baseEvent];
        plot_ciCohHistogram(sigciCohChanges_flatten, chnPairNames, f_selected, titlename, 'histClim', [-1 1],...
            'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1]);
        saveimgname = [saveimg_prefix '-' pdcond '-b' baseEvent  '-' compEvent '_align2' char(align2_comp) '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        close(gcf)
        
        %%% clear
        clear pdcond ciCohChangesfile
        clear ciCohChanges psedoiCohChanges f_selected ntrials
        clear sigciCohChanges sigciCohChanges_flatten chnPairNames
        clear titlename saveimgname nshuffle
        clear align2_comp tdur_base align2name_comp align2name_base align2_base
    end
    
    clear compEvent baseEvent     
end



function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_goodTrials_align2(files, align2, tdur_trial, varargin)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
%         Args:
%             align2: the event to be aligned (e.g Event.REACHONSET or uint 1-5)
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%   
%           Name-Value: 
%               'codesavefolder' - code saved folder
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


coli_align2 = uint32(align2);

coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);

t_minmax_reach = 0.2;

% load fs_lfp and T_chnsarea
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
        if t_reach < t_minmax_reach 
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
        
        clear t_reach t_return idxdur lfp_phase_1trial
    end
    
    clear filename file 
    clear('lfpdata', 'T_idxevent_lfp', 'goodTrials');
end



function [lfptrials, fs_lfp, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, varargin)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately
% [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2PeakV(files, [t_AOI(1) t_AOI(2)], 'codesavefolder', savecodefolder);
%
%   not include trials with t_reach <0.2s
% 
%         Args:
%             align2: the event to be aligned 
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

coli_align2 = uint32(align2);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);

t_minmax_reach = 0.2;

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'selectedTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    ntrials = length(lfpdata);
    for tri = 1: ntrials
        
        % ignore trials marked with 0
        if ~selectedTrials(tri)
            continue
        end
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_minmax_reach 
            clear t_reach
            continue
        end
        
        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
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
        
        
        % extract phase for 1 trial
        lfp_1trial = lfpdata{tri};
        idxdur = round(tdur_trial * fs_lfp) + idx_time0;
        if idxdur(1) == 0
            idxdur(1) = 1;
        else
            idxdur(1) = idxdur(1) + 1;
        end
        lfp_phase_1trial = lfp_1trial(:,idxdur(1) :idxdur(2));
           
        % cat into lfptrials
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach idxdur lfp_phase_1trial lfp_1trial
    end
end





