function m3_uNHP_SKT_PSD(animal, varargin)
%%   PSD estimates for mild and normal individually
%
%       psd for each brain area, as well as each DBS contact
%
%

%%

% parse params
p = inputParser;
addParameter(p, 'newsavefolder', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
newsavefolder = p.Results.newsavefolder;

%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code') - 1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder, 'util')));
addpath(genpath(fullfile(codefolder, 'NHPs')));

[~, codefilename]= fileparts(codefilepath);
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'SKT','fs500Hz', codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);

%%  input setup

% variables for plotting
plotF_AOI = [8 40];

% input folder: extracted raw rest data with grayMatter
if strcmpi(animal, 'Kitty')
    inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');
end
if strcmpi(animal, 'Jo')
    inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');
end

ylimits = [0 0.25];

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
if(exist(savefolder, 'dir') && newsavefolder)
    rmdir(savefolder, 's')
end
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename_prefix = 'psd_';

%% Code Start Here
EventPhases = {'preMove';'earlyReach'};

if(strcmpi(animal, 'Jo'))
    cond_cell = {'normal', 'PD'};
end
if(strcmpi(animal, 'Kitty'))
    cond_cell = {'normal', 'moderate'};
end

chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
t_min_reach = 0.5;
for ei = 1: length(EventPhases)
    eventname = EventPhases{ei};
    
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI4PSD_extract(eventname, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
        
        if strcmpi(animal, 'Jo') && strcmpi(pdcond, 'PD')
            files1 = dir(fullfile(inputfolder, ['*_mild_*.mat']));
            files2 = dir(fullfile(inputfolder, ['*_moderate_*.mat']));
            files = [files1; files2];
            clear files1 files2
        else
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
        end
        
        if strcmpi(animal, 'Kitty')
            [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, t_AOI, 't_min_reach', t_min_reach);
        end
        if strcmpi(animal, 'Jo')
            [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, t_AOI, 't_min_reach', t_min_reach); % lfptrials: nchns * ntemp * ntrials
        end
        
        % remove unused chns
        ChnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
        lfptrials = lfptrials(ChnsOfI_mask, :, :);
        T_chnsarea = T_chnsarea(ChnsOfI_mask, :);
        
        % psd calculate
        [pxx, f_selected] = pxx_eacharea(lfptrials, T_chnsarea, fs, plotF_AOI);
        
        if strcmpi(pdcond, 'moderate')
            pdcond = 'PD';
        end
        
        pxxs.(eventname).(pdcond) = pxx;
        
        if ~exist('brainareas', 'var')
            brainareas = fieldnames(pxx);
        end
        
        clear pdcond align2 t_AOI
        clear files lfptrials fs T_chnsarea
        clear ChnsOfI_mask lfptrials T_chnsarea
        clear pxx
    end
end
save(fullfile(savefolder, 'psd.mat'), 'pxxs', 'f_selected', 'brainareas');


%%%  plot  %%%
for ei = 1: length(EventPhases)
    eventname = EventPhases{ei};
    for i = 1:length(brainareas)
        brainarea = brainareas{i};
        
        % load normal and PD data
        psd_normal = pxxs.(eventname).normal.(brainarea);
        psd_PD = pxxs.(eventname).PD.(brainarea);
        
        %  psd_normal, psd_mild: nfs * nsegs
        plotPSD_comp_1chn(psd_normal, psd_PD, f_selected, plotF_AOI, savefolder, brainarea, animal, eventname)
        
        clear brainarea psd_normal psd_PD
    end
end

close all



function [lfptrials, fs_lfp, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, tdur_trial, varargin)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
% Args:
%       files: the files used for extraction       
%       align2: the event to be aligned 
%       tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])     
%
%             
%
%       Name-Value: 
%            't_min_reach' - the minimal reach time used, default 0.5 s
%            't_max_reach' - the max reach time used, default inf
% 
% Return:
%             lfptrials: nchns * ntemp * ntrials
%             chnAreas:
%             fs: sampling rate


% parse params
p = inputParser;
addParameter(p, 't_min_reach', 0.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 't_max_reach', inf, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
t_min_reach = p.Results.t_min_reach;
t_max_reach = p.Results.t_max_reach;


% code start here
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
        if t_reach < t_min_reach || t_reach > t_max_reach
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
            
            if t_reachonset2peakV < 0.3 || t_peakV2reach < 0.3
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


function [lfptrials, fs_lfp, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, tdur_trial, varargin)
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
%
%       Name-Value: 
%            't_min_reach' - the minimal reach time used, default 0.5 s
%
%   
%           't_max_reach' - the max reach time used, default inf
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


% parse params
p = inputParser;
addParameter(p, 't_min_reach', 0.5, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 't_max_reach', inf, @(x) assert(isnumeric(x) && isscalar(x)));


parse(p,varargin{:});
t_min_reach = p.Results.t_min_reach;
t_max_reach = p.Results.t_max_reach;


% code start here
coli_align2 = uint32(align2);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);


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
        if t_reach < t_min_reach  || t_reach > t_max_reach
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
            
            if t_reachonset2peakV < 0.3 || t_peakV2reach < 0.3
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


function [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI4PSD_extract(event, animal, varargin)
%
%   
% 
%   Inputs:
%       animal
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
%           'pdcond' - if animal == Kitty, required;


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'pdcond', '', @isstr);
parse(p,varargin{:});

codesavefolder = p.Results.codesavefolder;
pdcond = p.Results.pdcond;

if strcmpi(animal, 'Kitty') && isempty(pdcond)
    disp('pdcond is required for animal Kitty')
    return;
end

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

if strcmpi(event, 'preMove')
    align2 = SKTEvent.TargetOnset;
    t_AOI = [-1 -0.8];
end

if strcmpi(event, 'earlyReach')
    align2 = SKTEvent.ReachOnset;
    t_AOI = [0 0.2];
end

align2name = char(align2);

% special case for Kitty
if strcmpi(animal, 'Kitty') && strcmpi(pdcond, 'normal') && strcmpi(event, 'preMove')
    align2 = SKTEvent.TargetOnset;
    t_AOI = [0.2 0.4];
    align2name = 'StartTrial';
end





function [psds, f_selected] = pxx_eacharea(lfptrials, T_chnsarea, fs, f_AOI)
%% extract psd of all the segments from lfptrials, each psd for each dbs contact and one psd for one area (except dbs)
%
% Inputs:
%       lfptrials: nchns * ntemp * ntrials

%
%
% Outputs:
%
%       psds: the PSD estimate of all the segments from the files
%           e.g. pxxs =
%                   struct with fields:
%                      M1: [nfs * nsegs double]
%                     stn0_1: [nfs * nsegs double]
%                      gp0_1: [nfs * nsegs double]
%
%       f_selected: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)



[nchns, ntemp, ntrials] = size(lfptrials);

psds = struct();

for chi = 1:nchns
    
    brainarea = T_chnsarea.brainarea{chi};
    
    % change stn0-1 to STN, gp0_1 to GP
    brainarea_saved = brainarea;
    if contains(brainarea_saved, 'stn')
        brainarea_saved = 'STN';
    end
    if contains(brainarea_saved, 'gp')
        brainarea_saved = 'GP';
    end
    for tri = 1:ntrials
        
        % extract the lfp data of segi , lfp_oneseg: 1* ntemp
        lfp_onetrial = squeeze(lfptrials(chi, :, tri));
        
        
        %%% calcualte PSD throuch zscore and pwelch %%%
        
        % zscore of the lfp_oneseg along each column
        x = zscore(lfp_onetrial);
        
        % fft
        [~, fx, ~, px] = spectrogram(x, ntemp, 0,[],fs);
        if  chi == 1 && tri == 1
            freqs = fx;
            idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
            f_selected = round(freqs(idx_f), 2);
            
            clear freqs
        end

        if ~isfield(psds, brainarea_saved) % first time calculate pxx for brainarea
            psds.(brainarea_saved) = px;
        else % f_selected: nfs * 1; psds.M1: nfs * nsegs
            psds.(brainarea_saved) = cat(2, psds.(brainarea_saved), px);
        end
        
        
        clear lfp_onetrial x fx ps
    end
    
    clear brainarea brainarea_saved
end

function plotPSD_comp_1chn(psd_normal, psd_PD, F_all, plotF_AOI, savefolder, brainarea, animal, eventname)
%%  plot the psd comparison of normal and mild
%
%   Inputs:
%       psd_normal, psd_mild: psd of all segments in normal or mild, nfs * nsegs
%
%       F_all: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
%
%       plotF_AOI: the plot frequence range (in hertz)  (1 * 2)
%
%       savefile_prefix: the savefile prefix including folder

% find the idx for F_AOI
idx_AOI = find(F_all >= plotF_AOI(1) & F_all <= plotF_AOI(2));

% colors setup
color_normal_range = [224, 255, 255] / 255;
color_normal_mean = [0, 0, 255] / 255;
color_PD_range = [238, 238, 238] / 255;
color_PD_mean = [0, 0, 0] / 255;

% plot setup
linewidth = 1.5;

% F_AOI
F_AOI = F_all(idx_AOI);
% reshape F_AOI
[m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m * n); clear m n

% extract the  psd of AOI
psd_normal_FAOI = psd_normal(idx_AOI, :);
psd_moderate_FAOI = psd_PD(idx_AOI, :);

psd_normal_high = max(psd_normal_FAOI, [], 2);
psd_normal_low = min(psd_normal_FAOI, [], 2);
psd_normal_mean = mean(psd_normal_FAOI, 2);

psd_PD_high = max(psd_moderate_FAOI, [], 2);
psd_PD_low = min(psd_moderate_FAOI, [], 2);
psd_PD_mean = mean(psd_moderate_FAOI, 2);

% reshape, psd_*_high/low/mean into 1 * nfs

[m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m * n); clear m n
[m, n] = size(psd_normal_low); psd_normal_low = reshape(psd_normal_low, 1, m * n); clear m n
[m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m * n); clear m n


[m, n] = size(psd_PD_high); psd_PD_high = reshape(psd_PD_high, 1, m * n); clear m n
[m, n] = size(psd_PD_low); psd_PD_low = reshape(psd_PD_low, 1, m * n); clear m n
[m, n] = size(psd_PD_mean); psd_PD_mean = reshape(psd_PD_mean, 1, m * n); clear m n

% plot range
figure
set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
fill([F_AOI flip(F_AOI)], [psd_normal_high flip(psd_normal_low)], color_normal_range, 'LineStyle', 'none')
hold all
fill([F_AOI flip(F_AOI)], [psd_PD_high flip(psd_PD_low)], color_PD_range, 'LineStyle', 'none')

% plot mean
h1 = plot(F_AOI, psd_normal_mean, 'Color', color_normal_mean, 'LineWidth', linewidth);
h3 = plot(F_AOI, psd_PD_mean, 'Color', color_PD_mean, 'LineWidth', linewidth);


xlim([min(F_AOI) max(F_AOI)])
ylim([0 0.2])
set(gca,'XTick',[10 15 20 25 30 35 40 45 50],'YTick',[0 0.1 0.2]);


xlabel('Frequency (Hz)', 'FontWeight','bold')
ylabel('Power', 'FontWeight','bold')

% legend
legend([h1, h3], {'Normal',  'PD'})


% save figure
savename = fullfile(savefolder, [animal eventname '_SKT_psd_' brainarea]);
print(gcf, savename, '-painters', '-depsc')
print(gcf, savename, '-dpng', '-r1000')



clear psd_allsegs_normal  psd_allsegs_PD
clear psd_normal_FAOI  psd_PD_FAOI
clear psd_normal_high psd_normal_low psd_normal_mean psd_PD_high psd_PD_low psd_PD_mean
clear h1 h2 h3 maxPSD F_maxPSD idx_max
clear savename

