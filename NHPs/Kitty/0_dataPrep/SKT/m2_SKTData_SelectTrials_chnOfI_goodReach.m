function m2_segSKTData_SelectTrials_chnOfI_goodReach(varargin)
%  plot spectrogram of each trial and manually select trials only use good
%  reach trials
%
%   Input:
%       Name-Value: 
%           animal
%           fi_str - select file start index, default = 1
%           fi_end - select file end index, default = length(files)
%           autoSaveMode - auto save mode, 'y' or 'n'(default)

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));



% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
global animal
animal = animal_extract(codecorresfolder);

%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_avgArea');

cond_cell = cond_cell_extract(animal);
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);


coli_targetonset = uint32(SKTEvent.TargetOnset);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);
coli_returnonset = uint32(SKTEvent.ReturnOnset);
coli_mouth = uint32(SKTEvent.Mouth);


imFormat = 'tif';

% align to event
align2 = SKTEvent.ReachOnset;
coli_align2 = uint32(align2);

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

savefolder_trials = fullfile(savefolder, 'trials');
if ~exist(savefolder_trials, 'dir')
    mkdir(savefolder_trials);
end

savefolder_aveDay = fullfile(savefolder, 'avgDay');
if ~exist(savefolder_aveDay, 'dir')
    mkdir(savefolder_aveDay);
end

%% Code Start Here
warning off
close all

files = dir(fullfile(inputfolder, '*.mat'));

% parse params
p = inputParser;
addParameter(p, 'fi_str', 1, @isscalar);
addParameter(p, 'fi_end', length(files), @isscalar);
addParameter(p, 'autoSaveMode', 'n', @(x) isscalar(x)&&ischar(x));
parse(p,varargin{:});
fi_str = p.Results.fi_str;
fi_end = p.Results.fi_end;
autoSaveMode =  p.Results.autoSaveMode;
if ~strcmpi(autoSaveMode, 'y')
    autoSaveMode =  'n';
end



[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', savecodefolder);


cond_cell = cond_cell_extract(animal);
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);


filei = fi_str;
while(filei <=  fi_end)
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'fs_lfp', 'mask_goodreach', 'T_chnsarea', 'T_idxevent_lfp', ...
        'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    
    ntrials = size(T_idxevent_lfp, 1);
    if ntrials == 1
        
        clear('lfpdata', 'fs_lfp', 'mask_goodreach', 'T_chnsarea', 'T_idxevent_lfp', ...
            'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
        filei = filei + 1;
        continue;
    end
    
    %%% zscore the lfp data
    for tri = 1 : length(lfpdata)
        lfp_1trial = lfpdata{tri};
        zscored_lfpdata = zeros(size(lfp_1trial));
        for chi = 1: size(lfp_1trial, 1)
            tmp = squeeze(lfp_1trial(chi, :));
            zscored_lfpdata(chi, :) = zscore(tmp);
            clear tmp
        end
        lfpdata{tri} = zscored_lfpdata;
        clear lfp_1trial zscored_lfpdata
    end
   
    
    % extract dateofexp, bktdt and pdcond
    idx_cond = cellfun(@(x) contains(lower(filename), x), cond_cell);
    pdcond = cond_cell{idx_cond};
    datestrmatch = regexp(filename, '_[0-9]*_', 'match');
    dateofexp = datenum(datestrmatch{1}(2:end-1), 'yyyymmdd');
    bkstrmatch = regexp(filename, '_bktdt[0-9]*', 'match');
    if isempty(bkstrmatch)
        nstrmatch = regexp(filename, '_n[0-9]', 'match');
        bkstr = ['n' nstrmatch{1}(3:end)];
        clear nstrmatch
    else
        bktdt = str2num(bkstrmatch{1}(7:end));
        bkstr = ['bktdt' num2str(bktdt)];
        clear bktdt
    end
    
    ntrials = length(lfpdata);
    savefile_selectedTrialsMarkers = fullfile(savefolder, [animal '_chnOfI_TrialsWMarkers_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_' bkstr  '.mat']);
    
    
    %%%--- show spectrogram of each trial and select goodTrials ----%%%%
    disp(['checking fi = '  num2str(filei) '/' num2str(fi_end) ' ' filename])
    
    % extract data of Chns of interest
    mask_usedChns = cellfun(@(x) contains(x, chnsOfI), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';
    for tri = 1 : length(lfpdata)
        lfpdata{tri} = lfpdata{tri}(mask_usedChns, :);
    end
    
    
    % reload saved markers if existed
    if exist(savefile_selectedTrialsMarkers, 'file')
        load(savefile_selectedTrialsMarkers, 'selectedTrials');
    else
        selectedTrials = true(ntrials, 1);
    end

    showname = ['fi = ' num2str(filei) ':' animal '-' pdcond '-' datestr(dateofexp, 'yyyymmdd') '-' bkstr];
    
    madata = smoothWspeed_trial;
    madata2 = Wrist_smooth_trial;
    maName = 'Wspeed  XYZ-wrist';
    
    % check spectrogram trial by trial
    selectedTrials = check_chnOfI_spectrogram(filename, lfpdata, T_idxevent_lfp, T_chnsarea, fs_lfp,... 
        madata, madata2, maName, fs_ma, T_idxevent_ma,showname, ...
        'mask_goodReach', mask_goodreach, 'inSelectedTrials', selectedTrials, 'savefolder', savefolder_trials, ...
        'fig_left', 200, 'fig_bottom', 200,  'fig_width', 1200, 'fig_height', 600, ...
        'autoSaveMode', autoSaveMode);


    save(savefile_selectedTrialsMarkers, 'selectedTrials', 'lfpdata', 'T_idxevent_lfp', 'fs_lfp', 'T_chnsarea',...
        'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    

    
    %%% --- show averaged one day spectrogram using selectedTrials--- %%%
    if strcmp(pdcond, 'normal')
        tdur_trial = tdur_trial_normal;
        t_minmax_reach = t_minmax_reach_normal;
    end
    if strcmp(pdcond, 'mild')
        tdur_trial = tdur_trial_mild;
        t_minmax_reach = t_minmax_reach_mild;
    end
    if strcmp(pdcond, 'moderate')
        tdur_trial = tdur_trial_moderate;
        t_minmax_reach = t_minmax_reach_moderate;
    end
    
    disp(['only average across trials reach time larger than ' num2str(t_minmax_reach(1)) ' s'])
    disp(['tdur =  [' num2str(tdur_trial(1)) ' ' num2str(tdur_trial(2)) '] s'])
    
    % extract phase trials respected to coli_align2
    lfp_phase_trials = [];
    for tri = 1: length(lfpdata)
        
        % ignore trials marked with 0
        if ~selectedTrials(tri)
            continue
        end
        
        % select trials based on reach and return duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_minmax_reach(1)
            disp(['trial i = ' num2str(tri) ', t_reach = ' num2str(t_reach)])
            clear t_reach
            continue
        end
        
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
        lfp_1trial = lfpdata{tri};
        lfp_phase_1trial = lfp_1trial(:, idxdur(1):idxdur(2)); % lfp_phase_1trial: nchns * ntemp
        
        lfp_phase_trials = cat(3, lfp_phase_trials, lfp_phase_1trial); % lfp_phase_trials: nchns * ntemp * ntrials
        
        clear t_reach t_return idxdur lfp_phase_1trial lfp_1trial
    end
    
    oneday_spectrogram_img = fullfile(savefolder_aveDay, [animal '_selectedTrials_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_' bkstr]);
    if isempty(lfp_phase_trials)
        disp('lfp_phase_trials is empty, skip spectrogram across trials')
        if exist(oneday_spectrogram_img, 'file') % delete if exist already
            delete(oneday_spectrogram_img);
        end
        filei = filei + 1;
        continue;
    end
    
    % plot spectrogram across all good trials of one day
    disp(['# of selected Trials for averaged spectrogram across trials: ' num2str(size(lfp_phase_trials, 3))])
    plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea,  tdur_trial, fs_lfp, animal, pdcond, align2, showname)
    clear idxGroupNames
    
    if ~strcmpi(autoSaveMode, 'y')
        % Recheck today or Check the next day
        reply = input(['Check the next day (y) or Recheck this day (n)[y]:'], 's');
        if isempty(reply) || lower(reply) ~= 'n'
            
            saveas(gcf, oneday_spectrogram_img, imFormat);
            
            filei = filei + 1;
        end
        
    else
        saveas(gcf, oneday_spectrogram_img, imFormat);
        filei = filei + 1;
    end
    
    
    close all
    
    
    clear filename lfpdata T_idxevent fs T_chnsarea tbl_goodTrialsMarks
    clear mask_STN mask_GP mask_Others idxGroups
    clear idx_cond datestrmatch pdcond dateofexp bktdt
    clear tdur_trial t_minmax_reach t_minmax_return
    clear lfpdata_origin T_idxevent_origin
end
end



function plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, tdur_trial, fs, animal, pdcond, align2, showname)
% plot lfpdata of all the channels: nchns * ntemp * ntrials

twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);

% subplot/Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;


subp_startLeft = 0.05;


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
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    psd_plot = psd(idx_f, :);
    times_plot = times + tdur_trial(1);
    
    % convert into dB and then gauss filted
    psd_plot = 10 * log10(psd_plot);
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[1,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot
end

% plot spectrogram
[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);
fig = figure('PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]); 
annotation(fig,'textbox',...
            [subp_startLeft 0.1  1 0.03],...
            'String', {showname}, ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
        
nchns = size(psd_allchns, 3);
for chi = 1 : nchns
    brainarea = T_chnsarea.brainarea{chi};
    
    % identify clim
    if contains(brainarea, 'stn')
        clim = clim_Spectrogram_STN;
    else
        if contains(brainarea, 'gp')
            clim = clim_Spectrogram_GP;
        else
            clim = clim_Spectrogram_Others;
        end
    end
    
    
    % plot
    subplot(nchns, 1, chi);
    imagesc(times_plot, freqs_plot, psd_allchns(:, :, chi)); hold on
    
    if ~exist('clim', 'var') || isempty(clim)
        set(gca,'YDir','normal')
    else
        set(gca,'YDir','normal', 'CLim', clim)
    end
    colormap(jet)
    colorbar
    
    % xlabel
    if chi == nchns
        xlabel('time/s')
        
        xtls = xticklabels;
        xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
        xticklabels(xtls)
        clear xtls
    else
        xticks([])
    end
    
    % title and ylabel
    ylabel('Frequency(Hz)')
    title([animal ' ' pdcond ':' brainarea])
    
    % plot reach onset line
    plot([0 0], ylim, 'r--', 'LineWidth',1.5)
end
end


function selectedTrials = check_chnOfI_spectrogram(filename, lfpdata, T_idxevent_lfp, T_chnsarea, fs_lfp, madata, madata2, maName, fs_ma, T_idxevent_ma, showname, varargin)
% 
%
%   Input:
%       lfpdata: nchns * ntemp * ntrials
%       T_idxevent_lfp
%       T_chnsarea
%       fs_lfp
%
%       Name-Value: 
%           'autoSaveMode' - auto save mode, 'y' or 'n'(default 'n')
%           'fig_left' - figure position left (default 50)
%           'fig_bottom' - figure position bottom (default 50)
%           'fig_width' - figure position width (default 1200)
%           'fig_height' - figure position height (default 60)
%           'inSelectedTrials' - input logical array for each trial selected (1) or not (0), 
%                                default set to be all ones
%           'mask_goodReach' -- input array for each trial good (1) or not (0), default ([])
%           'savefolder' -- savefolder, default pwd
%           't_showmaxbef' - max show duration before align to (negative, default [], i.e show all) 
%           't_showmaxaft' - max show duration after align to (positive, default [], i.e show all) 
%
%      Return:
%           selectedTrials: logical array for each trial selected (1) or
%           not (0), ntrials * 1 logical

ntrials = length(lfpdata);

% parse params
p = inputParser;
addParameter(p, 'autoSaveMode', 'n', @(x) isscalar(x)&&ischar(x));
addParameter(p, 'fig_left', 50, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_bottom', 50, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_width', 1200, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_height', 600, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'inSelectedTrials', [], @(x) assert(isempty(x)||(islogical(x)&&isvector(x) &&size(x,1)==ntrials)));
addParameter(p, 'mask_goodReach', [], @(x) assert(isempty(x)||(isnumeric(x)&&isvector(x) &&length(x)==ntrials)));
addParameter(p, 'savefolder', pwd, @(x) assert(ischar(x)));
addParameter(p, 't_showmaxbef', [], @(x) assert(isempty(x)||(isnumeric(x) && isscalar(x) && x < 0)));
addParameter(p, 't_showmaxaft', [], @(x) assert(isempty(x)||(isnumeric(x) && isscalar(x)) && x > 0));

parse(p,varargin{:});
autoSaveMode =  p.Results.autoSaveMode;
fig_left = p.Results.fig_left;
fig_bottom = p.Results.fig_bottom;
fig_width = p.Results.fig_width;
fig_height = p.Results.fig_height;
inSelectedTrials = p.Results.inSelectedTrials;
mask_goodReach = p.Results.mask_goodReach;
savefolder = p.Results.savefolder;
t_showmaxbef =  p.Results.t_showmaxbef;
t_showmaxaft =  p.Results.t_showmaxaft;

if ~strcmpi(autoSaveMode, 'y')
    autoSaveMode =  'n';
end
if ~isempty(inSelectedTrials)
    selectedTrials = inSelectedTrials;
else
    selectedTrials = true(ntrials,1);
end

imFormat = 'tif';

% global parameters
twin = 0.2;
toverlap = 0.15;
f_AOI = [5 100];


coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_mouth = uint32(SKTEvent.Mouth);


coli_align2 = coli_reachonset;
align2 = 'reach onset';

areaName_left = 0.003;

% spectrogram subplot parameters
subp_startLeft = 0.1;
subp_endLeft = 0.95;
subp_startTop = 0.95;
subp_width = 0.2;
subp_height = 0.15;
supb_deltaX = 0.01;
supb_deltaY = 0.05;

% trial number per subplot
ntrials_pSubplot = floor((subp_endLeft - subp_startLeft - subp_width)/(subp_width + supb_deltaX )) + 1;



btn_width = 50;
btn_height = 30;

% Next button pos
btnNext_pos_left = (subp_endLeft) * fig_width;
btnNext_pos_bottom = (0.05) * fig_height;

% Prev button pos
btnPrev_pos_left = btnNext_pos_left;
btnPrev_pos_bottom = btnNext_pos_bottom + btn_height + 10;

% Finish Button pos
btnFinish_pos_left = (subp_endLeft) * fig_width;
btnFinish_pos_bottom = btnPrev_pos_bottom + btn_height + 20;

% buttons for check/uncheck trials on current page pos
btnCheck_width = 100;
btnCheckPage_pos_left = btnNext_pos_left - btnCheck_width - 50;
btnCheckPage_pos_bottom = btnNext_pos_bottom;
btnUncheckPage_pos_left = btnCheckPage_pos_left - btnCheck_width - 10;
btnUncheckPage_pos_bottom = btnNext_pos_bottom;

% buttons for check/uncheck buttons all trials pos
btnCheckAll_pos_left = btnUncheckPage_pos_left - btnCheck_width - 30;
btnCheckAll_pos_bottom = btnNext_pos_bottom;
btnUncheckAll_pos_left = btnCheckAll_pos_left - btnCheck_width - 10;
btnUncheckAll_pos_bottom = btnNext_pos_bottom;

% checkbox for showing/not showing threshold line for MA 
cbxShowThre_pos_left = btnNext_pos_left - 70;
cbxShowThre_pos_bottom = btnNext_pos_bottom + 450;
cbxShowThre_pos_width = 150;
cbxShowThre_pos_height = 20;

eventline_colors = ['c', 'r', 'g', 'y', 'k'];


% trial number checkbox parameters
cb_width = 150;
cb_height = 50;


% ngs: total number of sub-figures
nSubplots = ceil(ntrials / ntrials_pSubplot);

% extract dateofexp, bktdt and pdcond
global animal
cond_cell = cond_cell_extract(animal);
idx_cond = cellfun(@(x) contains(lower(filename), x), cond_cell);
pdcond = cond_cell{idx_cond};
datestrmatch = regexp(filename, '_[0-9]*_', 'match');
dateofexp = datenum(datestrmatch{1}(2:end-1), 'yyyymmdd');
bkstrmatch = regexp(filename, '_bktdt[0-9]*', 'match');
if isempty(bkstrmatch)
    nstrmatch = regexp(filename, '_n[0-9]', 'match');
    bkstr = ['n' nstrmatch{1}(3:end)];
    clear nstrmatch
else
    bktdt = str2num(bkstrmatch{1}(7:end));
    bkstr = ['bktdt' num2str(bktdt)];
    clear bktdt
end

[clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);

fig = figure('PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height], 'PaperPositionMode', 'auto'); 

hcbs = zeros(ntrials, 1);

% checkedAllGs: labels for marking whether all sub-figures are checked
checkedAllGs = zeros(nSubplots, 1);


%%% the first sub-figure (e.g triali = 1: subp_ntrials)
gi = 0;
clf(fig);
% add btn_next button
c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [btnNext_pos_left btnNext_pos_bottom btn_width btn_height]);
c_next.Callback = @btn_nextTrials;
title(['gi = ' num2str(gi + 1) '/' num2str(nSubplots)])

tri_str = 1;
tri_end = ntrials_pSubplot;
plot_spectrogram()
checkedAllGs(gi + 1) = 1;



%%% autoSaveMode or not
if strcmpi(autoSaveMode, 'y') % autoSaveMode
    
    nextFig_isLast = func_showNextTrials();
    while(~nextFig_isLast)
        % show next fig
        nextFig_isLast = func_showNextTrials();
    end
    func_finish();
else
    uiwait(fig);
end

    function plot_spectrogram()
        % plot spectrogram in fig
        
        annotation(gcf,'textbox',...
            [0.05 0.1 1 0.03],...
            'String', {showname}, ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
        
        % add checked this page and unChecked this page buttons
        c_checkedThisFig = uicontrol(fig, 'Style','pushbutton','String','checked this page', 'Position', [btnCheckPage_pos_left btnCheckPage_pos_bottom btnCheck_width btn_height]);
        c_checkedThisFig.Callback = @btn_checkedThisFig;
        c_unCheckedThisFig = uicontrol(fig, 'Style','pushbutton','String','unchecked this page', 'Position', [btnUncheckPage_pos_left btnUncheckPage_pos_bottom btnCheck_width btn_height]);
        c_unCheckedThisFig.Callback = @btn_uncheckedThisFig;
        
        % add checked this page and unChecked this page buttons
        c_checkedAll = uicontrol(fig, 'Style','pushbutton','String','checked all', 'Position', [btnCheckAll_pos_left btnCheckAll_pos_bottom btnCheck_width btn_height]);
        c_checkedAll.Callback = @btn_checkedAll;
        c_unCheckedAll = uicontrol(fig, 'Style','pushbutton','String','unchecked all', 'Position', [btnUncheckAll_pos_left btnUncheckAll_pos_bottom btnCheck_width btn_height]);
        c_unCheckedAll.Callback = @btn_uncheckedAll;
        
        
        % add checkbox for showing/not showing ma threshold = 30 line
        cbx_showMAThreshold = uicontrol('Style','checkbox','Value', 0, 'Position', [cbxShowThre_pos_left cbxShowThre_pos_bottom cbxShowThre_pos_width cbxShowThre_pos_height], 'String', 'Show Threshold = 30');
        cbx_showMAThreshold.Callback = @checkbox_showMAThreshold;
        
        for tri = tri_str: tri_end
            lfp_1trial = lfpdata{tri};
            
            coli = mod(tri,ntrials_pSubplot);
            if coli == 0
                coli = ntrials_pSubplot;
            end
            
            % left postion for tri
            subp_left = (coli -1) * (subp_width + supb_deltaX )+ subp_startLeft;
            
            % plot trial number Name and Checkbox
            pos_cb_trialN_left = (subp_left + subp_width/2) * fig_width - cb_width/2;
            pos_cb_trialN_bottom = (subp_startTop) * fig_height;
            str_Trial = ['triali = ' num2str(tri) '/' num2str(ntrials)];
            if ~isempty(mask_goodReach) && ~mask_goodReach(tri)
                str_Trial = [str_Trial '(bad)'];
                selectedTrials(tri) = selectedTrials(tri) && mask_goodReach(tri);
            end
            hcbs(tri) = uicontrol('Style','checkbox','Value',  selectedTrials(tri),...
                'Position', [pos_cb_trialN_left pos_cb_trialN_bottom cb_width cb_height],'String', str_Trial);
            set(hcbs(tri),'Callback',{@box_value});
            clear str_Trial val
            
            
            nchns = size(lfp_1trial, 1);
            for chi = 1 : nchns
                x = lfp_1trial(chi, 1:T_idxevent_lfp{tri, coli_mouth});
                
                % spectrogram
                [~, freqs, times, psd] = spectrogram(x, round(twin * fs_lfp), round(toverlap * fs_lfp),[],fs_lfp); % psd: nf * nt
                
                % select freqs_plot and corresponding psd
                idx_f = (freqs >= f_AOI(1) & freqs <= f_AOI(2));
                freqs_plot =  freqs(idx_f);
                psd_plot = psd(idx_f, :);
                % convert to db and gaussfilt
                psd_plot = 10 * log10(psd_plot);
                psd_plot = imgaussfilt(psd_plot,'FilterSize',5);
                times_plot = times - T_idxevent_lfp{tri, coli_align2} / fs_lfp;
                
                
                % identify clim
                brainarea = T_chnsarea.brainarea{chi};
                if contains(brainarea, 'stn')
                    clim = clim_Spectrogram_STN;
                else
                    if contains(brainarea, 'gp')
                        clim = clim_Spectrogram_GP;
                    else
                        clim = clim_Spectrogram_Others;
                    end
                end
                
                % spectrogram subplot
                subp_bottom = subp_startTop - subp_height - chi * (subp_height + supb_deltaY);
                subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
                imagesc(times_plot, freqs_plot, psd_plot, 'Tag', [num2str(tri) '-' num2str(chi)]);
                if exist('clim', 'var')|| ~isempty(clim)
                    set(gca,'YDir','normal', 'CLim', clim)
                else
                    set(gca,'YDir','normal')
                end
                clear brainarea clim
                
                colormap(jet)
                colorbar
                if (chi == nchns)
                    xlabel('time/s')
                    xtks = xticks();
                    xtklabels = xticklabels();

                    
                    % add xtick == 0.5 and - 0.5
                    xl = xlim;
                    xt = -1;
                    if isempty(find(xtks == xt))
                        idx_smaller = find(xtks < xt);
                        if ~isempty(idx_smaller) && xl(1) <= xt
                            xtks = [xtks(1:idx_smaller(end)) xt xtks(idx_smaller(end)+1:end)];
                            xtklabels = [xtklabels(1:idx_smaller(end)); {num2str(xt)}; xtklabels(idx_smaller(end)+1:end)];
                        end
                        clear idx_smaller
                    end
                    xt = 1;
                    if isempty(find(xtks == xt))
                        idx_smaller = find(xtks < xt);
                        if ~isempty(idx_smaller) && xl(2) >= xt
                            xtks = [xtks(1:idx_smaller(end)) xt xtks(idx_smaller(end)+1:end)];
                            xtklabels = [xtklabels(1:idx_smaller(end)); {num2str(xt)}; xtklabels(idx_smaller(end)+1:end)];
                        end
                        clear idx_smaller
                    end
                    
                    set(gca,'xtick',xtks,'XTickLabel',xtklabels,'TickLabelInterpreter','latex');
                    
                    clear xtks xtklabels xl xt
                else
                    
                    xticks([])
                end
                if(coli == 1)
                    ylabel('Frequency(Hz)')
                else
                    yticks([])
                end
                hold on
                % plot event lines
                for eventi = 1 : width(T_idxevent_lfp)
                    t_event = (T_idxevent_lfp{tri, eventi} - T_idxevent_lfp{tri, coli_align2}) / fs_lfp;
                    plot([t_event t_event], ylim, [eventline_colors(eventi) '--'], 'LineWidth',1.5)
                    clear t_event
                end
                clear eventi
                
                % adjust to [t_maxbef t_maxaft] if too long 
                xlim_t = get(gca, 'xlim');
                if ~isempty(t_showmaxbef)
                    if xlim_t(1) < t_showmaxbef
                        xlim_t(1) = t_showmaxbef;
                    end
                end
                if ~isempty(t_showmaxaft)
                    if xlim_t(2) > t_showmaxaft
                        xlim_t(2) = t_showmaxaft;
                    end
                end
                set(gca, 'xlim', xlim_t)
                clear xlim _t
                
                % extract spectrogram width and height using the first trial/first chn
                spect_axis_standerd = findobj('Tag', [num2str(tri_str) '-' num2str(1)]).Parent;
                spect_pos_standerd = get(spect_axis_standerd, 'Position');
                subp_width_standerd = spect_pos_standerd(3);
                subp_height_standerd = spect_pos_standerd(4);
                set(findobj('Tag', [num2str(tri) '-' num2str(chi)]).Parent, 'Position', [subp_left, subp_bottom, subp_width_standerd, subp_height_standerd])
                clear spect_axis_standerd spect_pos_standerd subp_width_standerd subp_height_standerd
                
                
                
                
                
                if chi == 1
                    % plot MA data in the first row
                    madata_1trial = madata{tri};
                    madata2_1trial = madata2{tri};
                    ma_WSpeed = madata_1trial(1:T_idxevent_ma{tri, coli_mouth});
                    ma_W_X =  squeeze(madata2_1trial(1:T_idxevent_ma{tri, coli_mouth}, 1));
                    ma_W_Y =  squeeze(madata2_1trial(1:T_idxevent_ma{tri, coli_mouth}, 2));
                    ma_W_Z =  squeeze(madata2_1trial(1:T_idxevent_ma{tri, coli_mouth}, 3));
                    clear madata_1trial madata2_1trial
                    
                    spect_axis = findobj('Tag', [num2str(tri) '-' num2str(1)]).Parent;
                    spect_pos = get(spect_axis, 'Position');
                    if tri == tri_str
                        subp_width_ma = spect_pos(3);
                        subp_height_ma = spect_pos(4);
                    end
                    spect_xlim = get(spect_axis, 'XLim');
                    subp_left_ma = spect_pos(1);
                    subp_bottom_ma = subp_startTop - subp_height;
                    times_plot_ma = ([1: length(ma_WSpeed)] - T_idxevent_ma{tri, coli_align2} )/ fs_ma;
                    subplot('Position', [subp_left_ma, subp_bottom_ma, subp_width_ma, subp_height_ma])
                    
                    
                    
                    % plot ma data
                    plot(times_plot_ma, ma_WSpeed, 'b'); hold on
                    plot(times_plot_ma, ma_W_X, 'r', times_plot_ma, ma_W_Y, 'g', times_plot_ma, ma_W_Z, 'k');
                    set(gca, 'XLim', spect_xlim);
                    xticks([])
                    clear  ma_WSpeed ma_W_X ma_W_Y ma_W_Z
                    
                    % plot MA threshold = 30
                    plot(xlim, [30 30], 'c--', 'Visible', 'off', 'Tag', ['MAThreshold-tri= ' num2str(tri)]);
                    
                    % plot event lines
                    for eventi = 1 : width(T_idxevent_ma)
                        t_event = (T_idxevent_ma{tri, eventi} - T_idxevent_ma{tri, coli_align2}) / fs_ma;
                        plot([t_event t_event], ylim, [eventline_colors(eventi) '--'], 'LineWidth',1.5)
                        clear t_event
                    end
                    clear eventi
                    
                    
                    if(coli == 1)
                        hl = legend({'WSpeed','W-X', 'W-Y', 'W-Z'},'AutoUpdate','off', 'Location', 'Best');
                        set(hl,'Position',[0.92 0.86 0.047 0.064],...
                            'AutoUpdate','off');
                        
                        clear h1
                    end
                    
                    annotation(gcf,'textbox',...
                        [areaName_left subp_bottom_ma + subp_height / 2 0.03 0.03],...
                        'String', maName, ...
                        'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
                    
                    clear x_ma spect_pos spect_xlim subp_bottom_ma subp_left_ma  times_plot_ma
                end
                
                
                % plot chn Names only in the first shown trial
                if(coli == 1)
                    annotation(gcf,'textbox',...
                        [areaName_left subp_bottom + subp_height / 2 0.03 0.03],...
                        'String',T_chnsarea{chi, 'brainarea'}, ...
                        'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
                end
                
                
                clear x freqs times psd idx_f freqs_plot psd_plot times_plot
                clear subp_bottom
            end
            
            clear lfp_1trial
        end
        save_oneSubplotFig();
        
    end

    function btn_checkedThisFig(src, event)
        for tri = tri_str: tri_end
            set(hcbs(tri),'Value',1);
             selectedTrials(tri) = 1;
        end
    end

    function btn_uncheckedThisFig(src, event)
        for tri = tri_str: tri_end
            set(hcbs(tri),'Value',0);
            selectedTrials(tri) = 0;
        end
        
    end

    function btn_checkedAll(src, event)
        for tri = 1: ntrials
            selectedTrials(tri) = 1;
        end
        
        for tri = tri_str: tri_end
            set(hcbs(tri),'Value',1);
        end
    end

    function btn_uncheckedAll(~, event)
      
        for tri = 1: ntrials
            selectedTrials(tri) = 0;
        end
        
        for tri = tri_str: tri_end
            set(hcbs(tri),'Value',0);
        end
    end

    function save_oneSubplotFig()
        trials_spectrogram_img = fullfile(savefolder, [animal '_trials_spect_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_' bkstr '_trial' num2str(tri_str) '-' num2str(tri_end)]);
        saveas(gcf, trials_spectrogram_img, imFormat);
    end

    function nextFig_isLast = showNewTrials()
        clf(fig);
        
        if(gi < nSubplots-1)
            c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [btnNext_pos_left btnNext_pos_bottom btn_width btn_height]);
            c_next.Callback = @btn_nextTrials;
        end
        
        if(gi > 0)
            c_prev = uicontrol(fig, 'Style','pushbutton', 'String', 'Previous', 'Position', [btnPrev_pos_left btnPrev_pos_bottom btn_width btn_height]);
            c_prev.Callback = @btn_prevTrials;
        end
        
        % start and end trial number for next graph
        tri_str = gi * ntrials_pSubplot + 1;
        tri_end = (gi + 1) * ntrials_pSubplot;
        nextFig_isLast = false;
        if tri_end >= ntrials
            tri_end = ntrials;
            nextFig_isLast = true;
        end
        plot_spectrogram()
        
        checkedAllGs(gi + 1) = 1;
        if all(checkedAllGs)
            c_Finish = uicontrol(fig, 'Style','pushbutton', 'String', 'Finish', 'Position', [btnFinish_pos_left btnFinish_pos_bottom btn_width btn_height]);
            c_Finish.Callback = @btn_finish;
        end
        
    end


    function nextFig_isLast = func_showNextTrials()
        % save current figure and show next trials
        
        gi = mod(gi + 1, nSubplots);
        
        save_oneSubplotFig();
        
        nextFig_isLast = showNewTrials();
    end
    
    function func_showPrevTrials()
        % save current figure and show previous trials
        
        gi = mod(gi - 1, nSubplots);
        
        save_oneSubplotFig();
        
        showNewTrials();
    end

    function func_finish()
        % save current figure and then close it
       
        save_oneSubplotFig();
        
        close(fig)
    end

    function btn_finish(src,event)
        func_finish();
    end % finishCheck

    function btn_nextTrials(src,event)
        func_showNextTrials();
    end 

    function btn_prevTrials(src,event)
        func_showPrevTrials();
    end

    function box_value(hObj,event)
        % Called when boxes are used
        
        % find the used box
        idx = find(hcbs==hObj);
        
        selectedTrials(idx) =  get(hObj,'Value');
    end

    function checkbox_showMAThreshold(hObj, event)
        showMAThreshold = get(hObj, 'Value');
        

        if showMAThreshold 
            for tri = tri_str : tri_end
                set(findobj('Tag', ['MAThreshold-tri= ' num2str(tri)]), 'Visible', 'on');
            end
        else
            for tri = tri_str : tri_end
                set(findobj('Tag', ['MAThreshold-tri= ' num2str(tri)]), 'Visible', 'off');
            end
        end
    end

end
