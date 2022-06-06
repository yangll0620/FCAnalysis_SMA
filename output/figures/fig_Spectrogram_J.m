function fig_Spectrogram_J()
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
input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');
input_folder_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_segSKTData_SelectTrials_chnOfI');

input_restfile_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'Rest', 'm4_restData_PSD', 'psd__allsegs_normalmildmoderate.mat');
input_restfile_K = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'Rest', 'm2_restData_PSD', 'psd__allsegs_normalmildmoderate.mat');
chnsused_K = {'M1', 'stn1_2', 'gp1_2'};
chnsused_J = {'M1', 'stn1_2', 'gp3_4'};
plotF_AOI = [8 40];

tdur_trial_normal_J = [-0.8 0.8];
tdur_trial_mild_J = [-0.8 0.8];
tdur_trial_moderate_J = [-0.8 0.8];

t_min_reach_K = 0.2;
tdur_trial_normal_K = [-0.6 1];
tdur_trial_moderate_K = [-0.6 1];
t_min_reach_J = 0.2;
tdur_trial_normal_J = [-0.6 1];
tdur_trial_mild_J = [-0.6 1];
tdur_trial_moderate_J = [-0.6 1];

clims_J.STN = [-30 -15];
clims_J.GP = [-30 -15];
clims_J.M1 = [-30 -10];
clims_K.STN =  [-40 0];
clims_K.GP =  [-40 0];
clims_K.M1 =  [-40 0];
    


savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;


% align to event
align2 = SKTEvent.ReachOnset;
coli_align2 = uint32(align2);


%% plot figure parameters
w_colormap = 240; % width  for the colormap
h_colormap = 120; % height for the colormap

w_delta_RestSK = 40;
w_deltax_J = 30; % x distance between two color map for animal J
w_deltax_K = 20; % x distance between two color map for animal K

h_deltay_J = 10;
h_deltay_K = 15;

w_textChname = 80; % width showing the channel name, i.e. M1
w_textPowerLabel = 40; % height showing the frequency label, i.e. Power
w_textPowerNum = 10; % height showing the frequency number, i.e. 0, 0.1, 0.2
w_textFreLabel = 40; % height showing the frequency label, i.e. Frequences(Hz)
w_textFreNum = 10; % height showing the frequency number, i.e. 10 12
w_textColorbar = 80; % width showing the colarbar 


h_textCond = 30; % height showing the condition, i.e. Mild-Normal
h_textTimeNum = 10; % height showing the time number, i.e. -0.5 0 0.5
h_textTimeLabel = 40; % height showing the time label, i.e. time(s)
h_textFreNum = 10; % height showing the Freq number, i.e. 10 12
h_textFreLabel = 40; % height showing the frequency label, i.e. Frequences(Hz)

%% Code start here
cond_cell_J = cond_cell_extract('Jo');
cond_cell_K = cond_cell_extract('Kitty');


animal = 'Jo';
cond_cell = cond_cell_J;
inputfolder = input_folder_J;
input_restfile = input_restfile_J;
tdur_trial_normal = tdur_trial_normal_J;
tdur_trial_mild = tdur_trial_mild_J;
tdur_trial_moderate = tdur_trial_moderate_J;
t_min_reach = t_min_reach_J;
w_deltax = w_deltax_J;
h_deltay = h_deltay_J;
clims = clims_J;
chnsused = chnsused_J;

nrows = 3;
ncols = length(cond_cell);

close all
fig_width = w_textChname + (w_textPowerLabel + w_textPowerNum + w_colormap + w_delta_RestSK) +(w_textFreLabel + w_textFreNum) + w_colormap * ncols + w_deltax * (ncols-1) + w_textColorbar; 
fig_height = h_textCond  + h_colormap * nrows + h_deltay * (nrows-1) + (h_textTimeLabel + h_textTimeNum); 
fig = figure('Position', [50 50 fig_width fig_height]);
set(fig, 'PaperUnits', 'points');

%%% plot Rest PSD in column 1
load(input_restfile, 'pxxs_allfiles_normal', 'pxxs_allfiles_moderate', 'F_pxx');
idx_AOI = find(F_pxx >= plotF_AOI(1) & F_pxx <= plotF_AOI(2));
freqs = F_pxx(idx_AOI);

%pxxs_allfiles_moderate = pxxs_allfiles_moderate(idx_AOI, :);
nrows = length(chnsused);
w_left_start = w_textChname;
for rowi = 1 : nrows
    chnname = chnsused{rowi};
    
    psds.normal = pxxs_allfiles_normal.(chnname);
    psds.normal = psds.normal(idx_AOI, :);
    if exist('pxxs_allfiles_mild', 'var')
        psds.mild = pxxs_allfiles_mild.chnname;
        psds.mild = psds.mild(idx_AOI, :);
    end
    psds.moderate = pxxs_allfiles_moderate.M1;
    psds.moderate = psds.moderate(idx_AOI, :);
    
    show_FreqNum = false;
    show_FreqLabel = false;
    if rowi == nrows
        show_FreqNum = true;
        show_FreqLabel = true;
    end

    h_outer_top = h_textCond + (h_colormap + h_deltay) * (rowi -1);
    h_outer_diff = h_colormap;
    h_inner_top = 0;
    h_inner_bottom = 0;
    if show_FreqNum
        h_outer_diff = h_outer_diff + h_textFreNum;
        h_inner_bottom = h_inner_bottom + h_textFreNum;
    end
    if show_FreqLabel
        h_outer_diff = h_outer_diff + h_textFreLabel;
        h_inner_bottom = h_inner_bottom + h_textFreLabel;
    end
    h_outer_bottom = fig_height - (h_outer_top + h_outer_diff);
    
    w_outer_left = w_left_start;
    w_inner_left = w_textPowerLabel + w_textPowerNum;
    w_inner_right = 0;
    w_outer_diff = w_colormap + w_textPowerLabel + w_textPowerNum;
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);
    
    % outer and inner margin
    subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
    subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
    
    
    plotPSD_comp_1chn(psds, freqs, ...
    'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin, ...
    'show_xlabel', show_FreqLabel, 'show_xticklabels', show_FreqNum);

    clear chnname psds 
    clear show_FreqNum show_FreqLabel
    clear w_outer_left w_outer_right w_inner_left w_inner_right w_outer_diff 
    clear h_outer_top h_outer_bottom h_outer_diff h_inner_top h_inner_bottom
    clear subplot_outerMargin subplot_innerposMargin
end




%%% plot Spectrogram on the right
w_left_start = w_textChname + + (w_textPowerLabel + w_textPowerNum + w_colormap + w_delta_RestSK);
for coli = 1 : ncols
    pdcond = cond_cell{coli};
    
    %%% extract lfptrials
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));   
    eval(['tdur_trial = tdur_trial_' pdcond ';']);
    if strcmpi(animal, 'Kitty')
        [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, t_min_reach);
    end
    if strcmpi(animal, 'Jo')      
        [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, tdur_trial, t_min_reach);
    end
       
    %%% calc psd_allchns
    [psd_allchns, freqs, times] = calc_spectrogram(lfptrials, fs, tdur_trial);
    
    %%% set show inf and position
    show_freqLabel = false;
    show_freqNum = false;
    show_colorbar = false;
    if coli == 1
        show_freqLabel = true;
        show_freqNum = true;
    end
    if coli == ncols
        show_colorbar = true;
    end
    
    w_outer_left = w_left_start + (w_textFreLabel + w_textFreNum) + (w_colormap + w_deltax) * (coli -1);
    w_inner_left = 0;
    w_inner_right = 0;
    w_outer_diff = w_colormap;
    if show_freqLabel
        w_outer_left = w_outer_left - w_textFreLabel;
        w_outer_diff = w_outer_diff + w_textFreLabel;
        w_inner_left = w_inner_left + w_textFreLabel; 
    end
    if show_freqNum
        w_outer_left = w_outer_left - w_textFreNum;
        w_outer_diff = w_outer_diff + w_textFreNum;
        w_inner_left = w_inner_left + w_textFreNum; 
    end
    if show_colorbar
        w_outer_diff = w_outer_diff + w_textColorbar;
        w_inner_right = w_inner_right + w_textColorbar;
    end
    w_outer_right = fig_width - (w_outer_left + w_outer_diff);
    
    
    %%% plot
    for rowi = 1: nrows
        psd_1chn = squeeze(psd_allchns(:, :, rowi));
        
        % clim
        brainarea = T_chnsarea.brainarea{rowi};
        if contains(brainarea, 'stn')
            brainarea = 'STN';
        end
        if contains(brainarea, 'gp')
            brainarea = 'GP';
        end
        clim = clims.(brainarea);
        clear brainarea
        
        
        %%% show inf and position
        show_timeLabel = false;
        show_timeNum = false;
        if rowi == nrows
            show_timeLabel = true;
            show_timeNum = true;
        end
        
        h_outer_top = h_textCond + (h_colormap + h_deltay) * (rowi -1);
        h_outer_diff = h_colormap;
        h_inner_top = 0;
        h_inner_bottom = 0;
        if show_timeLabel
            h_outer_diff = h_outer_diff + h_textTimeLabel;
            h_inner_bottom = h_inner_bottom + h_textTimeLabel;
        end
        if show_timeNum
            h_outer_diff = h_outer_diff + h_textTimeNum;
            h_inner_bottom = h_inner_bottom + h_textTimeNum;
        end
        h_outer_bottom = fig_height - (h_outer_top + h_outer_diff);
        
        % outer and inner margin
        subplot_outerMargin = [w_outer_left h_outer_top w_outer_right h_outer_bottom];
        subplot_innerposMargin = [w_inner_left h_inner_top w_inner_right h_inner_bottom];
        
        
        % plot
        plot_spectrogram_1chn(psd_1chn, freqs, times, align2, clim, ...
            'fig', fig, 'outerposMargin', subplot_outerMargin, 'innerposMargin', subplot_innerposMargin,...
            'show_xlabel', show_timeLabel, 'show_xticklabels', show_timeNum, 'show_ylabel', show_freqLabel, 'show_yticklabels', show_freqNum, 'show_colorbar', show_colorbar);
        
        % clear
        clear psd_1chn clim
        clear show_timeLabel show_timeNum
        clear h_outer_top h_outer_bottom h_outer_diff h_inner_top h_inner_bottom
        clear subplot_outerMargin  subplot_innerposMargin
    end
    
    
    %%% final clear
    clear pdcond files tdur_trial lfptrials fs T_chnsarea
    clear psd_allchns freqs times
    clear show_freqLabel show_freqNum show_colorbar
    clear w_outer_left w_outer_right w_inner_left w_inner_right w_outer_diff
end
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
%           'innerposMargin' - [left top right bottom]  left, top, right, bottom of distance between outer and inner position default [5 5 5 5]
%           'outerposMargin' - [left top right bottom] outer pos margin between outer and figure bounder default [5 5 5 5]
%           'show_xticklabels' - show (true, default) or not show (false) xticklabels
%           'show_yticklabels' - show (true, default) or not show (false) yticklabels



% parse params
p = inputParser;
addParameter(p, 'fig', [], @(x) assert(isa(x, 'matlab.ui.Figure')));
addParameter(p, 'innerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'outerposMargin', [5 5 5 5], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'show_xlabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_xticklabels', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_ylabel', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'show_yticklabels', true, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
fig = p.Results.fig;
innerposMargin = p.Results.innerposMargin;
outerposMargin = p.Results.outerposMargin; 
show_xlabel = p.Results.show_xlabel;
show_xticklabels = p.Results.show_xticklabels;
show_ylabel = p.Results.show_ylabel;
show_yticklabels = p.Results.show_yticklabels;


%%% plot

% colors setup
color_range.normal = [224, 255, 255] / 255;
color_range.moderate = [238, 238, 238] / 255;
color_mean.normal = [0, 0, 255] / 255;
color_mean.moderate = [0, 0, 0] / 255;



if isempty(fig)
    fig = figure;
end

% set axes position
fig_pos = fig.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);
outerpos = [outerposMargin(1) outerposMargin(4) fig_width-outerposMargin(1)-outerposMargin(3) fig_height-outerposMargin(2)-outerposMargin(4)];
innerpos = [outerposMargin(1)+ innerposMargin(1) outerposMargin(4)+innerposMargin(4) outerpos(3)-innerposMargin(1)-innerposMargin(3) outerpos(4)-innerposMargin(2)-innerposMargin(4)];
ax = axes(fig, 'Units', 'pixels');
set(gca, 'OuterPosition', outerpos, 'Position', innerpos)



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
legend(hlegends, pdconds)

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

end







