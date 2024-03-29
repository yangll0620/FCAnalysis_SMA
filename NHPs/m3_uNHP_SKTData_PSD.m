function m3_uNHP_SKTData_PSD(animal, varargin)
%  extract lfp data respect to reachonset
% 
%  return:
%        lfptrials: nchns * ntemp * ntrials


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

% saved fig format
savefig_format = 'tif';

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



cond_cell = cond_cell_extract(animal);
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] ...
    = goodSKTTrials_reachReturn_tcritiria(animal);


tdur_trial = [-1.5 1];

% align to event
align2 = SKTEvent.ReachOnset;


f_AOI = [8 40];
t_AOI = [-0.5 0.5];
twin = 0.2;
toverlap = 0.18;





%% Code Start Here
EventPhases = {'preMove';'earlyReach'};

if(strcmpi(animal, 'Jo'))
    cond_cell = {'normal', 'PD'};
end
if(strcmpi(animal, 'Kitty'))
    cond_cell = {'normal', 'moderate'};
end

reply = input('save trials or no y/n [n] ', 's');
if isempty(reply) || ~strcmpi(reply, 'y')
    savetrials =  'n';
else
    savetrials = 'y';
end
clear reply

savefolder_trials = fullfile(savefolder, 'trials');
if ~exist(savefolder_trials, 'dir')
    mkdir(savefolder_trials);
end

chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
t_min_reach = 0.5;
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    
    eval(['t_minmax_reach = t_minmax_reach_' pdcond ';']);
    eval(['t_minmax_return = t_minmax_return_' pdcond ';']);
    
    % lfptrials
    [align2, t_AOI, ~] = SKT_EventPhase_align2_tAOI4PSD_extract(eventname, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
    if strcmpi(animal, 'Kitty')
        [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_K(files, align2, t_AOI, 't_min_reach', t_min_reach);
    end
    if strcmpi(animal, 'Jo')
        [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2_J(files, align2, t_AOI, 't_min_reach', t_min_reach); % lfptrials: nchns * ntemp * ntrials
    end
    
    
    % chnsofI
    mask_ChnsOfI = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_ChnsOfI, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    lfptrials = lfptrials(mask_ChnsOfI, :, :);
    
    %%% plot spectrogram across all trials  
    plot_spectrogram_acrossTrials(lfptrials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, ...
                                  'f_AOI', f_AOI, 't_AOI', t_AOI, 'twin', twin, 'toverlap', toverlap, ... 
                                  'savefolder', savefolder, 'savefig_format', savefig_format);
    saveas(gcf, fullfile(savefolder, [animal '_' pdcond]), savefig_format);
    close gcf
    
    if strcmpi(savetrials, 'y')
        %%% plot spectrogram for each trial
        for fi = 1 : length(files)
            file = fullfile(files(fi).folder, files(fi).name);
            plot_spectrogram_eachTrials(file, tdur_trial, align2, noisy_chns, savefolder_trials, savefig_format);
            clear file
        end
        
    end

    

    %%% final clear
    clear cond files tdur_trial t_minmax_reach t_minmax_return
    clear lfptrials fs T_chnsarea
   
end
end


function plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, varargin)
% plot lfpdata of all the channels: nchns * ntemp * ntrials


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [-0.5 0.5], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'twin', 0.2, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'toverlap', 0.18, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'savefolder', pwd, @isstr);
addParameter(p, 'savefig_format', '.tif', @isstr);


% parse parameters
parse(p,varargin{:});
f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
twin = p.Results.twin;
toverlap = p.Results.toverlap;
savefolder = p.Results.savefolder;
savefig_format = p.Results.savefig_format;



% subplot/Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 1800;
fig_height = 900;


subp_startLeft = 0.05;
subp_endLeft = 0.97;
subp_startTop = 0.98;
subp_width = 0.25;
subp_height = 0.12;
supb_deltaX = 0.02;
supb_deltaY = 0.015;

[psd_allchns, freqs_plot, times_plot] = calc_spectrogram_acrossTrials(lfp_phase_trials, tdur_trial, fs, ...
                                        'f_AOI', f_AOI, 't_AOI', t_AOI, 'twin', twin, 'toverlap', toverlap);

% plot spectrogram
% Group chns into STN, GP and others
mask_STN = contains(T_chnsarea.brainarea, 'stn');
mask_GP = contains(T_chnsarea.brainarea, 'gp');
mask_Others = ~(mask_STN | mask_GP);
idxGroups = [{find(mask_STN)}; {find(mask_GP)}; {find(mask_Others)}];
idxGroupNames = {'STN'; 'GP'; 'Others'};
[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal, 'codesavefolder', fullfile(savefolder, 'code'));
ntrials = size(lfp_phase_trials, 3);

fig = figure(); 
set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
for idxGi = 1 : length(idxGroups)
    idxs = idxGroups{idxGi};
    
    if strcmpi(idxGroupNames{idxGi}, 'STN')
        clim = clim_Spectrogram_STN;
    end
    if strcmpi(idxGroupNames{idxGi}, 'GP')
        clim = clim_Spectrogram_GP;
    end
    if strcmpi(idxGroupNames{idxGi}, 'Others')
        clim = clim_Spectrogram_Others;
    end
    
    
    subp_left = (idxGi -1) * (subp_width + supb_deltaX )+ subp_startLeft;
    
    for idxi = 1 : length(idxs)
        areai = idxs(idxi);
        brainarea = T_chnsarea.brainarea{areai};
        
        %%% plot subplot %%%
        subp_bottom = subp_startTop - subp_height - (idxi -1) * (subp_height + supb_deltaY);
        ax = axes(fig, 'Position', [subp_left, subp_bottom, subp_width, subp_height]);
        imagesc(ax, times_plot, freqs_plot, psd_allchns(:, :, idxs(idxi)));
        if isempty(clim)
            set(ax,'YDir','normal')
        else
            set(ax,'YDir','normal', 'CLim', clim)
        end
        
        colormap(jet)
        colorbar
        if idxi == length(idxs)
            xlabel('time/s')
  
            xtls = xticklabels(ax);
            xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
            xticklabels(xtls)
            clear xtls
        else
            xticks([])
        end
        
        ylabel('Frequency(Hz)')
        
        hold on
        % plot reach onset line
        plot(ax,[0 0], ylim, 'r--', 'LineWidth',1.5)
        
        title(ax, [animal ' ' pdcond ':' brainarea ', ntrials = ' num2str(ntrials)])
        clear subp_bottom ax 
        
        
        %%%  plot a sepately figure %%%
        fig_sep = figure(); set(fig_sep, 'PaperUnits', 'points',  'Position', [680 558 600 300]);
        ax_sep  = axes('Parent',fig_sep);
        imagesc(ax_sep, times_plot, freqs_plot, psd_allchns(:, :, idxs(idxi)));
        colormap(jet)
        colorbar('FontSize', 9)
        ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
        xlabel('time/s', 'FontSize', 112, 'FontWeight', 'bold')
        xtls = xticklabels(ax_sep);
        xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
        xticklabels(ax_sep, xtls)
        set(ax_sep,'fontsize',11)
        hold on
        % plot reach onset line
        plot(ax_sep, [0 0], ylim, 'r--', 'LineWidth',1.5)
        if isempty(clim)
            set(ax_sep,'YDir','normal')
        else
            set(ax_sep,'YDir','normal', 'CLim', clim)
        end
        set(ax_sep, 'Position', [0.09 0.15 0.8 0.79])
        title(ax_sep, [animal ' ' pdcond ':' brainarea ', ntrials = ' num2str(ntrials)]);
        savefile_sep = fullfile(savefolder, [animal '_sep_' pdcond '_' brainarea]);
        saveas(fig_sep, savefile_sep, savefig_format);
        close(fig_sep)
        
        clear fig_sep ax_sep  savefile_sep  xtls
        
        
        
        clear areai brainarea
    end
    
    clear clim
    
end
end


function [psd_allchns, freqs_plot, times_plot] = calc_spectrogram_acrossTrials(lfp_phase_trials, tdur_trial, fs, varargin)
% calc spectrogram of lfp_phase_trials: nchns * ntemp * ntrials
%   Inputs:
%       lfp_phase_trials: nchns * ntemp * ntrials
%
%   Returns:
%       psd_allchns: nf * nt * nchns
%       freqs_plot: nf * 1
%       times_plot: nt * 1


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [-0.5 0.5], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
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
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end

end


function [psd_allchns_alltrials, freqs_plot, times_plot] = calc_spectrogram_perTrial(lfp_phase_trials, tdur_trial, fs, varargin)
% calc spectrogram of lfp_phase_trials: nchns * ntemp * ntrials
%   Inputs:
%       lfp_phase_trials: nchns * ntemp * ntrials
%
%   Returns:
%       psd_allchns: nf * nt  * ntrials * nchns
%       freqs_plot: nf * 1
%       times_plot: nt * 1


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [8 40], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 't_AOI', [-0.5 0.5], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
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
psd_allchns_alltrials = [];
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
    
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psds(idx_f, idx_t, :);
    
    
    psd_allchns_alltrials = cat(4, psd_allchns_alltrials, psd_plot); % psd_allchns: nf * nt * ntrials * nchns
    clear psds psd psd_plot idx_t idx_f
end

end


function plot_spectrogram_eachTrials(file, tdur_trial, align2, noisy_chns, savefolder, varargin)
%
%   file: abs path for file

if length(varargin) >= 1
    savefig_format = varargin{1};
else
    savefig_format = 'tif';
end

coli_align2 = uint32(align2);


twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];
t_AOI = [-0.5 0.5];

idxGroupNames = {'STN'; 'GP'; 'M1'};

% subplot/Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 1800;
fig_height = 900;

subp_startLeft = 0.05;
subp_endLeft = 0.97;
subp_startTop = 0.98;
subp_width = 0.25;
subp_height = 0.12;
supb_deltaX = 0.02;
supb_deltaY = 0.015;


%%% code start here

[~,name,ext] = fileparts(file);
filename = [name, ext];
clear name ext


% extract pdcond datestring, bktdtstring and animal
if contains(filename, 'normal')
    pdcond = 'normal';
end
if contains(filename, 'mild')
    pdcond = 'mild';
end
if contains(filename, 'moderate')
    pdcond = 'moderate';
end

tmp = regexp(filename, '\d{8}', 'match');
datestring = tmp{1};
tmp = regexp(filename, 'bktdt\d*', 'match');
bktdtstring = tmp{1};
clear tmp

tmp = strfind(filename, '_');
animal = filename(1:tmp(1)-1);


% clim 
[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);


load(file, 'lfpdata', 'T_idxevent_lfp', 'fs_lfp', 'goodTrials', 'T_chnsarea',...
    'smoothWspeed_trial', 'Wrist_smooth_trial', 'T_idxevent_ma', 'fs_ma');


% remove the noisy chns
mask_noisyChns = cellfun(@(x) contains(x, noisy_chns), T_chnsarea.brainarea);
mask_notDBS_notM1 = ~strcmp(T_chnsarea.brainarea, 'M1') & ~strcmp(T_chnsarea.electype, 'DBS');
mask_usedChns = ~(mask_noisyChns | mask_notDBS_notM1);
T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';


        
nwin = round(twin * fs_lfp);
noverlap = round(toverlap * fs_lfp);  


% extract idxGroups for indices of each group
mask_STN = contains(T_chnsarea.brainarea, 'stn');
mask_GP = contains(T_chnsarea.brainarea, 'gp');
mask_Others = ~(mask_STN | mask_GP);
idxGroups = [{find(mask_STN)}; {find(mask_GP)}; {find(mask_Others)}];


% plot and save each trial
ntrials = length(goodTrials);
for tri = 1: ntrials
    
    % ignore trials marked with 0
    if ~goodTrials(tri)
        continue
    end
    
    
    %%% --- extract MA data: ma_WSpeed,  ma_W_X/Y/Z--- %%% 
    % extract trial with t_dur for MA data
    idxdur = round(tdur_trial * fs_ma) + T_idxevent_ma{tri, coli_align2};
    ma_WSpeed = smoothWspeed_trial(idxdur(1):idxdur(2), tri);
    ma_W_X =  squeeze(Wrist_smooth_trial(idxdur(1):idxdur(2), 1, tri));
    ma_W_Y =  squeeze(Wrist_smooth_trial(idxdur(1):idxdur(2), 2, tri));
    ma_W_Z =  squeeze(Wrist_smooth_trial(idxdur(1):idxdur(2), 3, tri));
    times_plot_ma = [1: length(ma_WSpeed)]/fs_ma + tdur_trial(1);
    clear idxdur
    
    
    
    
    %%% --- extract psd of lfp data for each chn: psd_allchns %%%    
    % extract trial with t_dur for lfp data: lfp_phase_1trial (nchns * ntemp)
    idxdur = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
    lfp_phase_1trial = squeeze(lfpdata(:, idxdur(1):idxdur(2), tri));
    lfp_phase_1trial = lfp_phase_1trial(mask_usedChns, :, :);
    clear idxdur
    
    % extract psd_allchns_plot using spectrogram, convert into dB and gaussfilt
    psd_allchns = []; % psd_allchns: nf * nt * nchns
    for chi = 1 : size(lfp_phase_1trial, 1)
        x = lfp_phase_1trial(chi, :);
        [~, freqs, times, psd] = spectrogram(x, nwin, noverlap,[],fs_lfp); % ps: nf * nt
        
        % convert into dB
        psd = 10 * log10(psd);
        
        % select freqs and corresponding psd
        idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        freqs_plot =  freqs(idx_f);
        times = times + tdur_trial(1);
        idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
        times_plot_lfp = times(idx_t);
        
        psd_plot = psd(idx_f, :);
        
        % cat into psd_allchns
        psd_allchns = cat(3, psd_allchns, psd_plot);
        
        clear x freqs times psd idx_f psd_plot idx_t
    end
    psd_allchns_plot = zeros(size(psd_allchns));
    for chi = 1: size(psd_allchns, 3)
        psd_allchns_plot(:, :, chi) = imgaussfilt(squeeze(psd_allchns(:, :, chi)),'FilterSize',[3,11]);
    end
    clear lfp_phase_1trial chi psd_allchns
    
    
    
    
    %%% --- plot and save spectrogram and MA data at the same time ---%%%
    fig = figure(); 
    set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
    for idxGi = 1 : length(idxGroups)
        idxs = idxGroups{idxGi};
        
        if strcmpi(idxGroupNames{idxGi}, 'STN')
            clim = clim_Spectrogram_STN;
            grp_name = 'STN';
        end
        if strcmpi(idxGroupNames{idxGi}, 'GP')
            clim = clim_Spectrogram_GP;
            grp_name = 'GP';
        end
        if strcmpi(idxGroupNames{idxGi}, 'M1')
            clim = clim_Spectrogram_Others;
            grp_name = 'M1';
        end
        
        % subp_left position 
        subp_left = (idxGi -1) * (subp_width + supb_deltaX )+ subp_startLeft;
        
        % plot for each chn in this group
        for idxi = 1 : length(idxs)
            areai = idxs(idxi);
            brainarea = T_chnsarea.brainarea{areai};
            
            subp_bottom = subp_startTop - subp_height - idxi * (subp_height + supb_deltaY);
            subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
            imagesc(times_plot_lfp, freqs_plot, psd_allchns_plot(:, :, idxs(idxi)), 'Tag', [grp_name '-' num2str(idxi)]);
            if isempty(clim)
                set(gca,'YDir','normal')
            else
                set(gca,'YDir','normal', 'CLim', clim)
            end
            
            colormap(jet)
            colorbar
            if idxi == length(idxs)
                xlabel('time/s')
                
                xtls = xticklabels;
                xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
                xticklabels(xtls)
                clear xtls
            else
                xticks([])
            end
            if strcmp(pdcond, 'normal')
                ylabel('Frequency(Hz)')
            else
                yticks([])
            end
            
            
            % plot reach onset line
            hold on;
            plot([0 0], ylim, 'r--', 'LineWidth',1.5)
            
            % title for this channel
            title([animal ' ' pdcond ':' brainarea ', triali = ' num2str(tri) '/' num2str(ntrials)])
            
            clear areai brainarea subp_bottom 
        end
        clear idxi
        
        %--- plot MA data in the first row
        spect_axis = findobj('Tag', [grp_name '-' num2str(1)]).Parent;
        spect_pos = get(spect_axis, 'Position');
        subp_left_ma = spect_pos(1);
        subp_bottom_ma = subp_startTop - subp_height;
        subp_width_ma = spect_pos(3);
        subp_height_ma = spect_pos(4);
        subplot('Position', [subp_left_ma, subp_bottom_ma, subp_width_ma, subp_height_ma])
        
        plot(times_plot_ma, ma_WSpeed, 'b'); hold on
        plot(times_plot_ma, ma_W_X, 'r', times_plot_ma, ma_W_Y, 'g', times_plot_ma, ma_W_Z, 'k');
        set(gca, 'XLim', get(spect_axis, 'XLim'));
        xticks([])
        
        % plot event lines
        plot([0 0], ylim, ['r--'], 'LineWidth',1.5)
        % plot legend once
        if idxGi == 1
            hl = legend({'WSpeed','W-X', 'W-Y', 'W-Z'});
            set(hl,'Position',[0.86 0.86 0.047 0.064],'AutoUpdate','off', 'Location', 'Best');
            
            clear h1
        end
        
        
        clear idxs clim grp_name subp_left
        clear spect_axis spect_pos subp_left_ma subp_bottom_ma subp_width_ma sup_height_ma
    end
    clear idxGi
    
    %%% save part
    savefile = fullfile(savefolder, [animal  '_' pdcond '_' datestring '_' bktdtstring '_trial' num2str(tri)]);
    saveas(gcf, savefile, savefig_format);
    close all
    
    %%% final clear
    clear ma_WSpeed ma_W_X ma_W_Y ma_W_Z times_plot_ma
    clear fig psd_allchns_plot time_plot_lfp freqs_plot
    clear savefile
end
end