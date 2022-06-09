function fig3_Spectrogram()

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

%% Input & save
[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
[~, funcname, ~]= fileparts(codefilepath);

input_SKTfolders.J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');
input_SKTfolders.K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');

input_restfiles.J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'Rest', 'm4_restData_PSD', 'psd__allsegs_normalmildmoderate.mat');
input_restfiles.K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'Rest', 'm2_restData_PSD', 'psd__allsegs_normalmildmoderate.mat');

chnsuseds.J = {'M1', 'stn1_2', 'gp3_4'};
chnsuseds.K = {'M1', 'stn1_2', 'gp1_2'};

plotF_AOI = [8 40];

tdur_trials_J.normal = [-0.8 0.8];
tdur_trials_J.mild = [-0.8 0.8];
tdur_trials_J.moderate = [-0.8 0.8];

tdur_trials_K.normal = [-0.6 1];
tdur_trials_K.moderate = [-0.6 1];


clims_J.STN = [-30 -15];
clims_J.GP = [-30 -15];
clims_J.M1 = [-30 -10];
clims_K.STN =  [-40 0];
clims_K.GP =  [-40 0];
clims_K.M1 =  [-40 0];


t_min_reach = 0.2;

savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end


savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;


%% Code start here

animals = {'Jo', 'Kitty'};

%%% plot Rest PSD
if false
    for ai = 1 : length(animals)
        animal = animals{ai};
        input_restfile = input_restfiles.(animal(1));
        chnsused = chnsuseds.(animal(1));
        
        RestPSD(input_restfile, plotF_AOI, chnsused, savefolder, animal);
        
        clear animal input_restfile chnsused
    end
end



%%% plot spectrogram
for ai = 1 : length(animals) 
    animal = animals{ai};
    inputfolder = input_SKTfolders.(animal(1));
    cond_cell = cond_cell_extract(animal);
    chnsused = chnsuseds.(animal(1));
    if strcmp(animal, 'Jo')
        tdur_trials = tdur_trials_J;
        clims = clims_J;
    end
    if strcmp(animal, 'Kitty')
        tdur_trials = tdur_trials_K;
        clims = clims_K;
    end

    %plotAllSpectrograms(cond_cell, inputfolder, animal, chnsused, tdur_trials, SKTEvent.ReachOnset, t_min_reach, clims, savefolder);
    
    if strcmp(animal, 'Kitty')
        plotAllSpectrograms(cond_cell, inputfolder, animal, chnsused, tdur_trials, SKTEvent.PeakV, t_min_reach, clims, savefolder);
    end
end



function plotAllSpectrograms(cond_cell, inputfolder, animal, chnsused, tdur_trials, align2, t_min_reach, clims, savefolder)
%
%
nconds = length(cond_cell);
for ci = 1 : nconds
    pdcond = cond_cell{ci};
    
    %%% extract lfptrials
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));   
    tdur_trial = tdur_trials.(pdcond);
    if strcmpi(animal, 'Kitty')
        [lfptrials, fs, T_chnsarea] = lfptrials_K_selectedTrials_align2(files, align2, tdur_trial, t_min_reach);
    end
    if strcmpi(animal, 'Jo')      
        [lfptrials, fs, T_chnsarea] = lfptrials_J_goodTrials_align2(files, align2, tdur_trial, t_min_reach);
    end
       
    %%% calc psd_allchns
    [psd_allchns, freqs, times] = calc_spectrogram(lfptrials, fs, tdur_trial);
    
    %%% set show inf and position
    show_freqLabel = false;
    show_freqNum = false;
    show_colorbar = false;
    if ci == 1
        show_freqLabel = true;
        show_freqNum = true;
    end
    if ci == nconds
        show_colorbar = true;
    end
    
       
    %%% plot
    T_chnsarea.brainarea = cellfun(@(x) strrep(x, '-', '_'), T_chnsarea.brainarea, 'UniformOutput', false);
    mask_used = cellfun(@(x) contains(x, chnsused), T_chnsarea.brainarea);
    psd_allchns = psd_allchns(:, :, mask_used);
    T_chnsarea = T_chnsarea(mask_used, :);
    nchns = size(psd_allchns, 3);
    for chi = 1: nchns
        psd_1chn = squeeze(psd_allchns(:, :, chi));
        
        % clim
        brainarea = T_chnsarea.brainarea{chi};
        ba = brainarea;
        if contains(ba, 'stn')
            ba = 'STN';
        end
        if contains(ba, 'gp')
            ba = 'GP';
        end
        clim = clims.(ba);
       
        
        
        %%% show inf and position
        show_timeLabel = false;
        show_timeNum = false;
        if chi == nchns
            show_timeLabel = true;
            show_timeNum = true;
        end
        
        % plot
        fig = figure('Position', [150 150 400 200]);
        plot_spectrogram_1chn(psd_1chn, freqs, times, align2, clim, ...
            'fig', fig, ...
            'show_xlabel', show_timeLabel, 'show_xticklabels', show_timeNum, 'show_ylabel', show_freqLabel, 'show_yticklabels', show_freqNum, 'show_colorbar', show_colorbar);
        
        subsavefile = ['Spectrogram-' animal '-' pdcond '-' brainarea];
        if align2 == SKTEvent.PeakV
            subsavefile = ['Spectrogram-' animal '-peakV-' pdcond '-' brainarea];
        end
        print(fig, fullfile(savefolder, subsavefile), '-painters', '-depsc');
        print(fig, fullfile(savefolder, subsavefile), '-dpng', '-r1000')
        close(fig)
        
        
        % clear
        clear psd_1chn clim
        clear show_timeLabel show_timeNum
        clear brainarea
    end
    
    
    %%% final clear
    clear pdcond files tdur_trial lfptrials fs T_chnsarea
    clear psd_allchns freqs times
    clear show_freqLabel show_freqNum show_colorbar
    clear w_outer_left w_outer_right w_inner_left w_inner_right w_outer_diff
end

function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_K_selectedTrials_align2(files, align2, tdur_trial, t_min_reach)
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
        if t_reach < t_min_reach 
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

function [lfptrials, fs_lfp, T_chnsarea] = lfptrials_J_goodTrials_align2(files, align2, tdur_trial, t_min_reach, varargin)
% extract lfp data respect to targetonset, reachonset, reach and returnonset separately

% 
%         Args:
%             align2: the event to be aligned (e.g Event.REACHONSET or uint 1-5)
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%             t_minmax_reach, t_minmax_return : min and max reach/return (s) for selecting trials (e.g [0.5 1])
%   
%           Name-Value: 
%               'codesavefolder' - code saved folder
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

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    filename = files(filei).name;
    file = fullfile(files(filei).folder, filename);
    
    load(file, 'lfpdata', 'goodTrials', 'T_idxevent_lfp', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
    
    
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
        if t_reach < t_min_reach(1) 
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
        lfp_1trial = squeeze(lfpdata(:,:,tri));
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
    
    clear filename file 
end


function [psd_allchns, freqs, times] = calc_spectrogram(lfp_phase_trials, fs, tdur_trial)
%
% Inputs:
%    lfp_phase_trials: nchns * ntemp * ntrials
%
% Return:
%   psd_allchns: nf * nt * nchns
%   freqs: nf * 1
%   times: 1 * nt
  
twin = 0.2;
toverlap = 0.18;
f_AOI = [8 40];
t_AOI = [-0.5 0.5];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);

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
    freqs =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end


function RestPSD(input_restfile, plotF_AOI, chnsused, savefolder, animal)
load(input_restfile, 'pxxs_allfiles_normal', 'pxxs_allfiles_moderate', 'F_pxx');
vars = whos('-file',input_restfile, 'pxxs_allfiles_mild');
if ~isempty(vars)
    load(input_restfile, 'pxxs_allfiles_mild');
end
idx_AOI = find(F_pxx >= plotF_AOI(1) & F_pxx <= plotF_AOI(2));
freqs = F_pxx(idx_AOI);


%pxxs_allfiles_moderate = pxxs_allfiles_moderate(idx_AOI, :);
nrows = length(chnsused);
show_PowerLabel = true;
show_PowerNum = true;
for rowi = 1 : nrows
    chnname = chnsused{rowi};
    
    psds.normal = pxxs_allfiles_normal.(chnname);
    psds.normal = psds.normal(idx_AOI, :);
    if exist('pxxs_allfiles_mild', 'var')
        psds.mild = pxxs_allfiles_mild.(chnname);
        psds.mild = psds.mild(idx_AOI, :);
    end
    psds.moderate = pxxs_allfiles_moderate.M1;
    psds.moderate = psds.moderate(idx_AOI, :);
    
    show_FreqNum = false;
    show_FreqLabel = false;
    show_legend = false;
    
    if rowi == 1
        show_legend = true;
    end
    if rowi == nrows
        show_FreqNum = true;
        show_FreqLabel = true;
    end
    
    fig = figure('Position', [150 150 300 200]);
    plotPSD_comp_1chn(psds, freqs,...
        'fig', fig,...
        'show_xlabel', show_FreqLabel, 'show_xticklabels', show_FreqNum, ...
        'show_ylabel', show_PowerLabel, 'show_yticklabels', show_PowerNum, 'show_legend', show_legend);
    
    subsavefile = ['RestPSD-' animal '-' chnname];
    print(fig, fullfile(savefolder, subsavefile), '-painters', '-depsc');
    print(fig, fullfile(savefolder, subsavefile), '-dpng', '-r1000')
    close(fig)
    

    clear chnname psds 
    clear show_FreqNum show_FreqLabel
end


function plotPSD_comp_1chn(psds, freqs, varargin)
% Inputs
%
%       psds.normal
%       psds.mild
%       psds.moderate
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default []
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default []
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels



% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'innerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_legend', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
fig = p.Results.fig;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin; 
show_xlabel = p.Results.show_xlabel;
show_xticklabels = p.Results.show_xticklabels;
show_ylabel = p.Results.show_ylabel;
show_yticklabels = p.Results.show_yticklabels;
show_legend = p.Results.show_legend;


%%% plot

% colors setup
color_range.normal = [224, 255, 255] / 255;
color_range.mild = [255, 228, 225] / 255;
color_range.moderate = [238, 238, 238] / 255;
color_mean.normal = [0, 0, 255] / 255;
color_mean.mild = [255, 00, 0] / 255;
color_mean.moderate = [0, 0, 0] / 255;



if isempty(fig)
    fig = figure();
end

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
ax = axes(fig, 'Units', 'pixels');
if ~isempty(outerposMargin) && ~isempty(innerposMargin)
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
end

pdconds = fieldnames(psds);
hlegends = [];
[m, n] = size(freqs); freqs = reshape(freqs, 1, m * n); clear m n
for pdi = 1 : length(pdconds)
    pdcond = pdconds{pdi};
    psd = psds.(pdcond);
    
    psd_high = max(psd, [], 2);
    psd_low = min(psd, [], 2);
    
    % reshape, psd_high/low into 1 * nfs   
    [m, n] = size(psd_high); psd_high = reshape(psd_high, 1, m * n); clear m n
    [m, n] = size(psd_low); psd_low = reshape(psd_low, 1, m * n); clear m n
    
    
    % plot
    fill(ax, [freqs flip(freqs)], [psd_high flip(psd_low)], color_range.(pdcond), 'LineStyle', 'none')
    hold all

    clear pdcond psd psd_high psd_low
end
for pdi = 1 : length(pdconds)
    pdcond = pdconds{pdi};
    psd = psds.(pdcond);
    
    psd_mean = mean(psd,2);
   
    % reshape, psd_mean into 1 * nfs   
    [m, n] = size(psd_mean); psd_mean = reshape(psd_mean, 1, m * n); clear m n

    h = plot(ax, freqs, psd_mean, 'Color', color_mean.(pdcond), 'LineWidth', 1);
    hlegends = [hlegends, h];
    
    clear pdcond psd psd_mean h
end

xlim([min(freqs) max(freqs)])
ylim([0 0.2])



% set show inf
if show_xticklabels
    set(ax,'XTick',[10 15 20 25 30 35 40 45 50]);
else
    set(ax,'XTick',[]);
end

if show_yticklabels
    set(ax,'YTick',[0 0.1 0.2]);
else
    set(ax,'YTick',[]);
end

if show_xlabel
    xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
end

if show_ylabel
    ylabel('Power', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
end

if show_legend
    legend(hlegends, pdconds)
end


function plot_spectrogram_1chn(psd_plot, freqs_plot, times_plot, align2, clim, varargin)
%
%   Inputs:
%
%       Name-Value: 
%           'fig' - figure handle to show the image (default [] to create a new one)
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default []
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default []
%           'show_xlabel' - show (true, default) or not show (false) xlabel
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_ylabel' - show (true, default) or not show (false) ylabel
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels
%           'show_colorbar' - show (true, default) or not show (false) colorbar
%



% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'outerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'innerposMargin', [], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_colorbar', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
fig = p.Results.fig;
outerposMargin = p.Results.outerposMargin; 
innerposMargin = p.Results.innerposMargin;
show_xlabel = p.Results.show_xlabel;
show_xticklabels = p.Results.show_xticklabels;
show_ylabel = p.Results.show_ylabel;
show_yticklabels = p.Results.show_yticklabels;
show_colorbar = p.Results.show_colorbar;



if isempty(fig)
    fig = figure;
end

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
ax = axes(fig, 'Units', 'pixels');
if ~isempty(outerposMargin) && ~isempty(innerposMargin)
    outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
    innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
    set(gca, 'OuterPosition', outerpos, 'Position', innerpos)
end


% plot
imagesc(ax, times_plot, freqs_plot, psd_plot); hold on
set(ax,'YDir','normal', 'CLim', clim)
colormap(jet)
c = colorbar;
c.Visible = 'off';


% plot reach onset line
plot(ax, [0 0], ylim, 'r--', 'LineWidth',1.5)

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

if show_ylabel
  pos = get(ax, 'Position');
  ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
  set(ax, 'Position', pos);
end

if show_yticklabels
else
    pos = get(ax, 'Position');
    yticks([]);
    set(ax, 'Position', pos);
end

if show_colorbar
  c.Visible = 'on';
  c.Label.String = 'Power (dB/Hz)';
end

