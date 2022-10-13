function fig3_spectrogram(varargin)
%
% fig3_spectrogram('newsavefolder', true);

% parse params
p = inputParser;
addParameter(p, 'newsavefolder', false, @(x) assert(islogical(x) && isscalar(x)));


parse(p,varargin{:});
newsavefolder = p.Results.newsavefolder;

codefilepath = mfilename('fullpath');


% find the codefolder
tmp = regexp(codefilepath, '.*\code', 'match');
if length(tmp) ~= 1
    disp('can not find code path correctly.')
    return;
end
codefolder = tmp{1};
clear tmp

% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));


%% Input
[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
[~, funcname, ~]= fileparts(codefilepath);


animals = {'Jo'; 'Kitty'};
cond_cell = {'normal', 'PD'};
EventPhases = {'preMove';'earlyReach'};

twin_plot.preMove.Jo = [-1 -0.5];
twin_plot.preMove.Kitty = [0.2 0.7];
twin_plot.earlyReach = [-0.5 1];

tdur_calc.preMove.Jo = [-1 0];
tdur_calc.preMove.Kitty = [0.2 1.2];
tdur_calc.earlyReach = [-1 1.5];

clims.Jo.STN = [-30 -15];
clims.Jo.GP = [-30 -15];
clims.Jo.M1 = [-30 -10];
clims.Kitty.STN = [-30 -5];
clims.Kitty.GP = [-30 -15];
clims.Kitty.M1 = [-30 -10];

t_min_reach = 0.5;

f_AOI = [8 40];


pos_ifigs.earlyReach = [50 50 360 150];
pos_ifigs.preMove = [50 50 150 150];

%% save setup
savefolder = fullfile(outputfolder, 'results', 'figures', 'current', funcname);
savecodefolder = fullfile(savefolder, 'code');
if(exist(savefolder, 'dir') && newsavefolder)
    rmdir(savefolder, 's')
end
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
if ~exist(savefolder, "dir")
    mkdir(savefolder);
end
copyfile2folder(codefilepath, savecodefolder);


% aisave folder
aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', 'Current',funcname);
if(exist(aisavefolder, 'dir') && newsavefolder)
    rmdir(aisavefolder, 's')
end
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end



%% Code Start Here
psdfile = fullfile(savefolder, 'psd.mat');

%%% extract psdfile
if ~exist(psdfile, 'file')
    for ai = 1 : length(animals)
        animal = animals{ai};
        
        chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
        
        if strcmpi(animal, 'Kitty')
            inputfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');
        end
        if strcmpi(animal, 'Jo')
            inputfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');
        end
        
        for ei = 1: length(EventPhases)
            eventname = EventPhases{ei};
            if strcmpi(eventname, 'preMove')
                tdur_trial = tdur_calc.(eventname).(animal);
            else
                tdur_trial = tdur_calc.(eventname);
            end
            for ci = 1 : length(cond_cell)
                pdcond = cond_cell{ci};
                
                if strcmpi(pdcond, 'PD')
                    files1 = dir(fullfile(inputfolder, ['*_mild_*.mat']));
                    files2 = dir(fullfile(inputfolder, ['*_moderate_*.mat']));
                    files = [files1; files2];
                    clear files1 files2
                else
                    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
                end
                
                [align2, ~, ~] = SKT_EventPhase_align2_tAOI4PSD_extract(eventname, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
                if strcmpi(animal, 'Kitty')
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, tdur_trial, 't_min_reach', t_min_reach);
                end
                if strcmpi(animal, 'Jo')
                    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, tdur_trial, 't_min_reach', t_min_reach); % lfptrials: nchns * ntemp * ntrials
                end
                
                
                % chnsofI
                ChnsOfI_mask = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
                lfptrials = lfptrials(ChnsOfI_mask, :, :);
                T_chnsarea = T_chnsarea(ChnsOfI_mask, :);
                clear ChnsOfI_mask
                
                % pxx
                [psd_allchns, fxs, txs] = calc_spectrogram_acrossTrials(lfptrials, tdur_trial, fs, 't_AOI', [], 'f_AOI', f_AOI);
                
                psds.(animal).(eventname).(pdcond) = psd_allchns;
                freqs.(animal).(eventname).(pdcond) = fxs;
                ts.(animal).(eventname).(pdcond) = txs;
                aligns2.(animal).(eventname).(pdcond) = align2;
                
                clear pdcond align2 lfptrials fs 
                clear psd_allchns fxs txs
            end
            
            clear eventname t_AOI tdur_trial ci
        end
        
        clear animal chnsOfI inputfolder ei
    end
    brainareas = {'M1', 'STN', 'GP'};
    save(psdfile, 'psds', 'freqs', 'ts', 'aligns2', 'brainareas');
    clear psds freqs ts align2 brainareas
end


%%% plot
load(psdfile, 'psds', 'freqs', 'ts', 'aligns2', 'brainareas');
animals = fieldnames(psds);
for ai = 1 : length(animals)
    animal = animals{ai};
    EventPhases = fieldnames(psds.(animal));
    for ei = 1: length(EventPhases)
        eventname = EventPhases{ei};
        cond_cell = fieldnames(psds.(animal).(eventname));
        for ci = 1 : length(cond_cell)
            pdcond = cond_cell{ci};
            psd_allchn = psds.(animal).(eventname).(pdcond);
            freqs_plot = freqs.(animal).(eventname).(pdcond);
            times_plot = ts.(animal).(eventname).(pdcond);
            align2 = aligns2.(animal).(eventname).(pdcond);
            if strcmpi(eventname, 'preMove')
                t_plot = twin_plot.(eventname).(animal);
            else
                t_plot = twin_plot.(eventname);
            end
            idx_t = (times_plot >= t_plot(1) & times_plot <= t_plot(2));
            times_plot = times_plot(idx_t);
            for chi = 1 : size(psd_allchn, 3)
                psd_plot = squeeze(psd_allchn(:, idx_t, chi));
                brainarea = brainareas{chi};
                
                savefilename = ['spect_' animal '_' brainarea  '_' pdcond '_' eventname];
                
                show_xticklabels = false;
                show_xlabel = false;
                show_titlename = false;
                
                if chi ==  size(psd_allchn, 3)
                    show_xticklabels = true;
                    show_xlabel = true;
                end
                
                show_colorbar = false;
                show_yticklabels = false;
                show_yticks = true;
                if ci == length(cond_cell) && ei == length(EventPhases)
                    show_colorbar = true;
                end
                if ci == 1 && ei == 1
                    show_yticklabels = true;
                end
                
                clim = clims.(animal).(brainarea);
                
                
                % plot
                pos_ifig = pos_ifigs.(eventname);
                plot_spectrogram_1chn(psd_plot, freqs_plot, times_plot, align2, clim, ...
                    pos_ifig, savefolder, savefilename, ...
                    show_xticklabels, show_xlabel, show_titlename, show_yticklabels, show_yticks, show_colorbar)
            end
            
            clear pdcond psd_allchn
        end
        clear eventname cond_cell
    end
    clear animal EventPhases
    
    clear animal ei
    
    close all
end

%%% copy to aifolder
copyfile2aifolder(savefolder, aisavefolder, 'eps');



function plot_spectrogram_1chn(psd_plot, freqs_plot, times_plot, align2, clim, ...
    pos_ifig, savefolder, savefilename, ...
    show_xticklabels, show_xlabel, show_titlename, show_yticklabels, show_yticks, show_colorbar)


ifig = figure('Position', pos_ifig);
set(ifig, 'PaperUnits', 'points');
ax = axes(ifig, 'Units', 'pixels');

% plot
imagesc(ax, times_plot, freqs_plot, psd_plot); hold on
colormap(jet)
set(ax,'YDir','normal', 'CLim', clim)


% plot vertical line x == 0
plot(ax, [0 0], ylim, 'r--', 'LineWidth',1.5)

c = colorbar;
c.Visible = 'off';
c.Label.String = 'Power (dB/Hz)';
climits = c.Limits;
cticks = (ceil(climits(1)/5)*5) : 5 : (floor(climits(2)/5)*5);
c.Ticks = cticks;



%%% show inf
if show_xlabel
  pos = get(ax, 'Position');
  xlabel('time (s)', 'FontSize', 12, 'FontWeight', 'bold')
  set(ax, 'Position', pos)
end

if show_xticklabels
    pos = get(ax, 'Position');

    xtls = xticklabels(ax);
    xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
    xticklabels(ax, xtls)
    set(ax, 'Position', pos);
else
    pos = get(ax, 'Position');
    xticks([]);
    set(ax, 'Position', pos);
end

if show_titlename
    title(titlename, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontname)
end

if show_yticklabels
  pos = get(ax, 'Position');
  ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
  set(ax, 'Position', pos);
end

if show_yticks
else
    pos = get(ax, 'Position');
    yticks([]);
    set(ax, 'Position', pos);
end

if show_colorbar
  c.Visible = 'on';
end

% save 
savefile = fullfile(savefolder, savefilename);
print(ifig, savefile, '-painters', '-depsc')
print(ifig, savefile, '-dpng', '-r1000')



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


function [psd_allchns, fxs, txs] = calc_spectrogram_acrossTrials(lfp_phase_trials, tdur_trial, fs, varargin)
% calc spectrogram of lfp_phase_trials: nchns * ntemp * ntrials
%   Inputs:
%       lfp_phase_trials: nchns * ntemp * ntrials
%
%   Returns:
%       psd_allchns: nf * nt * nchns
%       fxs: nf * 1
%       txs: 1 *nt


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [], @(x) assert(isempty(x)||(isnumeric(x) && isvector(x) && length(x)==2)));
addParameter(p, 'twin', 0.2, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'toverlap', 0.18, @(x) assert(isnumeric(x) && isscalar(x)));


% parse parameters
parse(p,varargin{:});
f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
twin = p.Results.twin;
toverlap = p.Results.toverlap;


%% Code Start here
nwin = round(twin * fs);
noverlap = round(toverlap * fs);


% calculate psd for each chn across trials
psd_allchns = [];
for chi = 1 : size(lfp_phase_trials, 1)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfp_phase_trials, 3)
        x = lfp_phase_trials(chi, :, tri);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    fxs =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    if ~isempty(t_AOI)
        idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
        txs = times(idx_t);
        psd_plot = psd(idx_f, idx_t);
    else
        txs = times;
        psd_plot = psd(idx_f, :);
    end
    
    
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd idx_f freqs_plot times psd_plot
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


