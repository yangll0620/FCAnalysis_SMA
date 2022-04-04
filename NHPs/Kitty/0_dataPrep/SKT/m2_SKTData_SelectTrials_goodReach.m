function m2_segSKTData_SelectTrials_goodReach()
%  plot spectrogram of each trial and manually select trials only use good
%  reach trials
%

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
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '/', '[A-Za-z]*']);
elseif ispc
    % Code to run on Windows platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '\\', '[A-Za-z]*']);
else
    disp('Platform not supported')
end
animal = codecorresfolder(fi + length('NHPs') + 1:j);

%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_avgArea');

cond_cell = cond_cell_extract(animal);
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);

if strcmpi(animal, 'bug')
    
    tdur_trial_normal = [-0.6 1];
    tdur_trial_mild = [-0.6 1];
    tdur_trial_moderate = [-0.6 1];
end
if strcmpi(animal, 'jo')

    tdur_trial_normal = [-0.8 0.8];
    tdur_trial_mild = [-0.8 0.8];
    tdur_trial_moderate = [-0.8 0.8];
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    
    tdur_trial_normal = [-0.6 1];
    tdur_trial_moderate = [-0.6 1];
end

if strcmpi(animal, 'pinky')
    
    tdur_trial_normal = [-0.6 1];
    tdur_trial_mild = [-0.6 1];
    tdur_trial_moderate = [-0.6 1];
end


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

%% starting
warning off
savefolder_trials = fullfile(savefolder, 'trials');
if ~exist(savefolder_trials, 'dir')
    mkdir(savefolder_trials);
end

savefolder_aveDay = fullfile(savefolder, 'avgDay');
if ~exist(savefolder_aveDay, 'dir')
    mkdir(savefolder_aveDay);
end


close all
files = dir(fullfile(inputfolder, '*.mat'));


% input start file number or start with 1
reply = input('Input start/end file (e.g.[20 30] or 20 for start) number = ', 's');
if isempty(reply)
    filei = 1; 
    nfiles = length(files);
else
    % parse start file number
    tmp = regexp(reply, '[*\d*', 'match');
    if ~isempty(tmp)
        startStr = tmp{1};
        if contains(startStr, '[')
            filei = str2num(startStr(2:end));
        else
            filei = str2num(startStr);
        end
        clear startStr
    else
        filei = 1; 
    end
    clear tmp
    
    % parse end file number
    tmp = regexp(reply, '\d*]', 'match');
    if ~isempty(tmp)
        endStr = tmp{1};
        nfiles = str2num(endStr(1:end-1));
        clear endStr
    else
        nfiles = length(files);
    end
    clear tmp
end
clear reply


% auto save without click
reply = input('Auto save mode [y] or no [n] ', 's');
if isempty(reply) || ~strcmpi(reply, 'y')
    autoSaveMode =  'n';
else
    autoSaveMode = 'y';
end
clear reply

while(filei <=  nfiles)
    
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
    savefile_selectedTrialsMarkers = fullfile(savefolder, [animal '_TrialsWMarkers_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_' bkstr  '.mat']);
    
    
    %%%--- show spectrogram of each trial in STN, GP and others groups and select goodTrials ----%%%%
    disp(['checking fi = '  num2str(filei) '/' num2str(nfiles) ' ' filename])
    
    % Group chns into STN, GP and others
    mask_STN = contains(T_chnsarea.brainarea, 'stn');
    mask_GP = contains(T_chnsarea.brainarea, 'gp');
    mask_Others = ~(mask_STN | mask_GP);
    idxGroups = [{find(mask_STN)}; {find(mask_GP)}; {find(mask_Others)}];
    groupNames = {'STN', 'GP', 'Others'};
    % split if more than 7 chns in one group
    idxGroups_split = {};
    groupNames_split = {};
    for idxGi = 1 : length(idxGroups)
        idxs = idxGroups{idxGi};
        if length(idxs) > 7
            splitN = ceil(length(idxs) / 7);
            for si = 1 : splitN
                i_str = (si - 1) * 7 + 1;
                i_end = si * 7;
                if i_end > length(idxs)
                    i_end = length(idxs);
                end
                idxGroups_split = [idxGroups_split; idxs(i_str:i_end)];
                groupNames_split = [groupNames_split; [groupNames{idxGi} num2str(si)]];
                clear i_str i_end
            end
        else
            idxGroups_split = [idxGroups_split; idxs];
            groupNames_split = [groupNames_split; groupNames{idxGi}];
        end
    end
    
    % reload saved markers if existed
    selectedTrials_Loaded = false;
    if exist(savefile_selectedTrialsMarkers, 'file')
        load(savefile_selectedTrialsMarkers, 'tbl_selectedTrialsMarks', 'idxGroups');
        groupNames_load = tbl_selectedTrialsMarks.Properties.VariableNames;
        idxGroups_load = idxGroups; clear idxGroups
        gName_iGroup_Equal = true;
        if length(groupNames_load) == length(groupNames_split) && length(idxGroups_split) == length(groupNames_split) && length(groupNames_load) == length(idxGroups_load)
            for gi = 1 :  length(groupNames_split)
                if ~strcmp(groupNames_split{gi}, groupNames_load{gi}) || ~isequal(idxGroups_split{gi}, idxGroups_load{gi})
                    gName_iGroup_Equal = false;
                    break
                end
            end
        else
            gName_iGroup_Equal = false;
        end
        
        if gName_iGroup_Equal
            selectedTrials_Loaded = true;
        end
    end
    idxGroups = idxGroups_split;
    groupNames = groupNames_split;
    if selectedTrials_Loaded
        
        selectedTrials_allGs = tbl_selectedTrialsMarks{:, :};
    else
        selectedTrials_allGs = ones(ntrials, length(idxGroups));
    end
    clear idxGroups_split groupNames_split


    showname = ['fi = ' num2str(filei) ':' animal '-' pdcond '-' datestr(dateofexp, 'yyyymmdd') '-' bkstr];
    
    % check_spectrogram_oneGroup
    [clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);
    for idxGi = 1 : length(idxGroups)
        idxs = idxGroups{idxGi};
        lfpdata_1group = cell(length(lfpdata), 1);
        for tri = 1: length(lfpdata)
            lfp_1trial = lfpdata{tri};
            lfpdata_1group{tri} = lfp_1trial(idxs, :);
            clear lfp_1trial
        end
        clear tri
        T_chnsarea_1group = T_chnsarea(idxs, :);
        
        
        if all(contains(T_chnsarea_1group.brainarea, 'stn'))
            clim_Spectrogram = clim_Spectrogram_STN;
            groupname = 'STN';
        else
            if all(contains(T_chnsarea_1group.brainarea, 'gp'))
                clim_Spectrogram = clim_Spectrogram_GP;
                groupname = 'GP';
            else
                clim_Spectrogram = clim_Spectrogram_Others;
                
                if any(contains(T_chnsarea_1group.brainarea, 'M1')) % only show M1, ignore SMA, VPLo et al.
                    M1_mask = contains(T_chnsarea_1group.brainarea, 'M1');
                    groupname = 'M1';
                    for tri = 1 : length(lfpdata_1group)
                        lfp_1trial = lfpdata_1group{tri};
                        lfpdata_1group{tri} = lfp_1trial(M1_mask, :);
                        clear lfp_1trial
                    end
                    T_chnsarea_1group = T_chnsarea_1group(M1_mask, :);
                else
                    continue;
                end    
                
            end
        end
              
        selectedTrials_1Grp = selectedTrials_allGs(:, idxGi);
        madata = smoothWspeed_trial;
        madata2 = Wrist_smooth_trial;
        maName = 'Wspeed  XYZ-wrist';
        selectedTrials_allGs(:, idxGi) = check_spectrogram_oneGroup(lfpdata_1group, T_idxevent_lfp, T_chnsarea_1group, fs_lfp, selectedTrials_1Grp,... 
            madata, madata2, maName, fs_ma, T_idxevent_ma,...
            showname, clim_Spectrogram, ...
            savefolder_trials, animal, groupname, pdcond, datestr(dateofexp, 'yyyymmdd'), bkstr, autoSaveMode, mask_goodreach);
        
        
        clear idxs lfpdata_1group T_chnsarea_1group goodTrials_1Grp
        
        
        %%% save selectedTrials and tbl_selectedTrialsMarks
        % selectedTrials = selectedTrials_allGs(:, 1) & selectedTrials_allGs(:, 2) & selectedTrials_allGs(:, 3)
        selectedTrials = ones(size(selectedTrials_allGs, 1), 1);
        for gi = 1 : size(selectedTrials_allGs, 2)
            selectedTrials = selectedTrials & selectedTrials_allGs(:, gi);
        end
        tbl_selectedTrialsMarks = array2table(selectedTrials_allGs, 'VariableNames', groupNames);
        if exist(savefile_selectedTrialsMarkers, 'file')
            delete(savefile_selectedTrialsMarkers)
        end
        save(savefile_selectedTrialsMarkers, 'selectedTrials', 'tbl_selectedTrialsMarks', 'idxGroups', 'lfpdata', 'T_idxevent_lfp', 'fs_lfp', 'T_chnsarea',...
            'fs_ma', 'T_idxevent_ma', 'smoothWspeed_trial', 'Wpos_smooth_trial', 'Wrist_smooth_trial');
    end
    

    
    
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
    
    
    % extract phase trials respected to coli_align2
    lfp_phase_trials = [];
    for tri = 1: length(lfpdata)
        
        % ignore trials marked with 0
        if ~selectedTrials(tri)
            continue
        end
        
        % select trials based on reach and return duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_minmax_reach(1) || t_reach > t_minmax_reach(2)
            disp(['trial i = ' num2str(tri) ', t_reach = ' num2str(t_reach)])
            selectedTrials(tri) = 0;
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
    idxGroupNames = tbl_selectedTrialsMarks.Properties.VariableNames;
    plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, idxGroups, idxGroupNames, tdur_trial, fs_lfp, animal, pdcond, align2, showname)
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



function plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, idxGroups, idxGroupNames, tdur_trial, fs, animal, pdcond, align2, showname)
% plot lfpdata of all the channels: nchns * ntemp * ntrials

twin = 0.2;
toverlap = 0.15;
f_AOI = [8 40];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);

% subplot/Figure parameters
fig_left = 2561;
fig_bottom = -300;
fig_width = 1920;
fig_height = 960;


subp_startLeft = 0.05;
subp_endLeft = 0.98;
subp_startTop = 0.98;
subp_width = 0.25;
subp_height = 0.12;
supb_deltaX = 0.02;
supb_deltaY = 0.015;
nGrps = length(idxGroups);
if (nGrps - 1) * supb_deltaX + nGrps * subp_width +  subp_startLeft > subp_endLeft
    subp_width = round((subp_endLeft - subp_startLeft - (nGrps - 1) * supb_deltaX)/nGrps, 2);
end


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
fig = figure(); set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
annotation(gcf,'textbox',...
            [subp_startLeft 0.1  1 0.03],...
            'String', {showname}, ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
for idxGi = 1 : length(idxGroups)
    idxs = idxGroups{idxGi};
    
    if contains(idxGroupNames{idxGi}, 'STN')
        clim = clim_Spectrogram_STN;
    end
    if contains(idxGroupNames{idxGi}, 'GP')
        clim = clim_Spectrogram_GP;
    end
    if contains(idxGroupNames{idxGi}, 'Others')
        clim = clim_Spectrogram_Others;
        idx_M1 = -1;
        for i = 1: length(idxs)
            if strcmp(T_chnsarea.brainarea{idxs(i)}, 'M1')
                idx_M1 = idxs(i);
            end
        end
        if idx_M1 == -1
            continue;
        else
            idxs = idx_M1;
        end
    end
    
    

    
    subp_left = (idxGi -1) * (subp_width + supb_deltaX )+ subp_startLeft;
    
    for idxi = 1 : length(idxs)
        areai = idxs(idxi);
        brainarea = T_chnsarea.brainarea{areai};
        
        subp_bottom = subp_startTop - subp_height - (idxi -1) * (subp_height + supb_deltaY);
        subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
        imagesc(times_plot, freqs_plot, psd_allchns(:, :, idxs(idxi)));
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
        
        ylabel('Frequency(Hz)')
        
        hold on
        % plot reach onset line
        plot([0 0], ylim, 'r--', 'LineWidth',1.5)
        
        
        title([animal ' ' pdcond ':' brainarea])
    end
    
    clear clim
    
end
end


function selectedTrials = check_spectrogram_oneGroup(lfpdata, T_idxevent_lfp, T_chnsarea, fs_lfp, selectedTrials, madata, madata2, maName, fs_ma, T_idxevent_ma, showname, clim, varargin)
% lfpdata: nchns * ntemp * ntrials


% parse varargin
if length(varargin) >= 8
    mask_goodReach = varargin{8};
else
    mask_goodReach = [];
end

if length(varargin) >= 7
    autoSaveMode = varargin{7};
else
    autoSaveMode = 'n';
end

if length(varargin) >= 6
    bkstr = varargin{6};
else
    bkstr = '';
end
if length(varargin) >= 5
    dateofexp_yyyymmdd = varargin{5};
else
    dateofexp_yyyymmdd = '';
end
if length(varargin) >= 4
    pdcond = varargin{4};
else
    pdcond = '';
end
if length(varargin) >= 3
    groupname = varargin{3};
else
    groupname = '';
end
if length(varargin) >= 2
    animal = varargin{2};
else
    animal = '';
end
if length(varargin) >= 1
    savefolder = varargin{1};
else
    savefolder = '';
end


imFormat = 'tif';

% global parameters
twin = 0.2;
toverlap = 0.15;
f_AOI = [5 100];

coli_targetonset = uint32(SKTEvent.TargetOnset);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);
coli_returnonset = uint32(SKTEvent.ReturnOnset);
coli_mouth = uint32(SKTEvent.Mouth);
% align to event
coli_align2 = coli_reachonset;
align2 = 'reach onset';

% subplot/Figure parameters
fig_left = 2561;
fig_bottom = -300;
fig_width = 1920;
fig_height = 960;


areaName_left = 0.003;

% spectrogram subplot parameters
subp_startLeft = 0.05;
subp_endLeft = 0.97;
subp_startTop = 0.95;
subp_width = 0.2;
subp_height = 0.12;
supb_deltaX = 0.01;
supb_deltaY = 0.01;
subp_nchns = floor((subp_startTop - subp_height)/(subp_height + supb_deltaY )) + 1;
subp_ntrials = floor((subp_endLeft - subp_startLeft - subp_width)/(subp_width + supb_deltaX )) + 1;


% position of Next and previous button
btn_width = 50;
btn_height = 30;
pos_btn_next_left = (subp_endLeft) * fig_width;
pos_btn_next_bottom = (0.05) * fig_height;
pos_btn_prev_left = pos_btn_next_left;
pos_btn_prev_bottom = pos_btn_next_bottom + btn_height + 10;
% Finish Button
pos_btn_finish_left = (subp_endLeft) * fig_width;
pos_btn_finish_bottom = 0.5 * fig_height;

eventline_colors = ['c', 'r', 'g', 'y', 'k'];


% trial number checkbox parameters
cb_width = 150;
cb_height = 50;


ntrials = length(lfpdata);

% ngs: total number of sub-figures 
ngs = ceil(ntrials / subp_ntrials);

fig = figure(); set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height], 'PaperPositionMode', 'auto');
hcbs = zeros(ntrials, 1);

% checkedAllGs: labels for marking whether all sub-figures are checked
checkedAllGs = zeros(ngs, 1);


%%% the first sub-figure (e.g triali = 1: subp_ntrials)
gi = 0;
clf(fig);
% add btn_next button
c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_btn_next_left pos_btn_next_bottom btn_width btn_height]);
c_next.Callback = @btn_nextTrials;
title(['gi = ' num2str(gi + 1) '/' num2str(ngs)])

tri_str = 1;
tri_end = subp_ntrials;
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

    function nextFig_isLast = showNewTrials()
        clf(fig);
        
        if(gi < ngs-1)
            c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_btn_next_left pos_btn_next_bottom btn_width btn_height]);
            c_next.Callback = @btn_nextTrials;
        end
        
        if(gi > 0)
            c_prev = uicontrol(fig, 'Style','pushbutton', 'String', 'Previous', 'Position', [pos_btn_prev_left pos_btn_prev_bottom btn_width btn_height]);
            c_prev.Callback = @btn_prevTrials;
        end
        
        % start and end trial number for next graph
        tri_str = gi * subp_ntrials + 1;
        tri_end = (gi + 1) * subp_ntrials;
        nextFig_isLast = false;
        if tri_end >= ntrials
            tri_end = ntrials;
            nextFig_isLast = true;
        end
        plot_spectrogram()
        
        checkedAllGs(gi + 1) = 1;
        if all(checkedAllGs)
            c_Finish = uicontrol(fig, 'Style','pushbutton', 'String', 'Finish', 'Position', [pos_btn_finish_left pos_btn_finish_bottom btn_width btn_height]);
            c_Finish.Callback = @btn_finish;
        end
        
    end


    function nextFig_isLast = func_showNextTrials()
        % save current figure and show next trials
        
        gi = mod(gi + 1, ngs);
        
        trials_spectrogram_img = fullfile(savefolder, [animal '_' groupname '_trials_spect_' pdcond '_' dateofexp_yyyymmdd '_' bkstr '_trial' num2str(tri_str) '-' num2str(tri_end)]);
        saveas(gcf, trials_spectrogram_img, imFormat);
        
        nextFig_isLast = showNewTrials();
    end
    
    function func_showPrevTrials()
        % save current figure and show previous trials
        
        gi = mod(gi - 1, ngs);
        
        trials_spectrogram_img = fullfile(savefolder, [animal '_' groupname '_trials_spect_' pdcond '_' dateofexp_yyyymmdd '_' bkstr '_trial' num2str(tri_str) '-' num2str(tri_end)]);
        saveas(gcf, trials_spectrogram_img, imFormat);
        
        showNewTrials();
    end

    function func_finish()
        % save current figure and then close it
       
        trials_spectrogram_img = fullfile(savefolder, [animal '_' groupname '_trials_spect_' pdcond '_' dateofexp_yyyymmdd '_' bkstr '_trial' num2str(tri_str) '-' num2str(tri_end)]);
        saveas(gcf, trials_spectrogram_img, imFormat);
        
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

    function plot_spectrogram()
        % plot lfpdata_1group of all the channels: nchns * ntemp * ntrial
       
        annotation(gcf,'textbox',...
            [0.87 0.45 1 0.03],...
            'String', {showname}, ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
        
        % add checked this page and unChecked this page buttons
        c_checkedThisFig = uicontrol(fig, 'Style','pushbutton','String','checked this page', 'Position', [1800 900 100 20]);
        c_checkedThisFig.Callback = @btn_checkedThisFig;
        c_unCheckedThisFig = uicontrol(fig, 'Style','pushbutton','String','unchecked this page', 'Position', [1800 870 100 20]);
        c_unCheckedThisFig.Callback = @btn_uncheckedThisFig;
        
        % add checked this page and unChecked this page buttons
        c_checkedAll = uicontrol(fig, 'Style','pushbutton','String','checked all', 'Position', [1800 800 100 20]);
        c_checkedAll.Callback = @btn_checkedAll;
        c_unCheckedAll = uicontrol(fig, 'Style','pushbutton','String','unchecked all', 'Position', [1800 770 100 20]);
        c_unCheckedAll.Callback = @btn_uncheckedAll;
        
        
        % add checkbox for showing/not showing ma threshold = 30 line
        cbx_showMAThreshold = uicontrol('Style','checkbox','Value', 0, 'Position', [1800 400 100 20], 'String', 'Show Threshold = 30');
        cbx_showMAThreshold.Callback = @checkbox_showMAThreshold;
        
        t_maxbef = -5;
        t_maxaft = 5;
        for tri = tri_str: tri_end
            lfp_1trial = lfpdata{tri};
            
            coli = mod(tri,subp_ntrials);
            if coli == 0
                coli = subp_ntrials;
            end
            
            % left postion for tri
            subp_left = (coli -1) * (subp_width + supb_deltaX )+ subp_startLeft;
            
            % plot trial number Name and Checkbox
            pos_cb_trialN_left = (subp_left + subp_width/2) * fig_width;
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
                
                                
                % spectrogram subplot
                subp_bottom = subp_startTop - subp_height - chi * (subp_height + supb_deltaY);
                subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
                imagesc(times_plot, freqs_plot, psd_plot, 'Tag', [num2str(tri) '-' num2str(chi)]);
                if ~isempty(clim)
                    set(gca,'YDir','normal', 'CLim', clim)
                else
                    set(gca,'YDir','normal')
                end
                
                colormap(jet)
                colorbar
                if (chi == nchns)
                    xlabel('time/s')
                    xtks = xticks();
                    xtklabels = xticklabels();
                    tmp = strsplit(align2);
                    xtklabels{xtks == 0} = ['$$\begin{array}{c}' ...
                               tmp{1} '\\'...
                               tmp{2} '\\'...
                               '\end{array}$$'];
                                        
                    % add xtick == 0.5 and - 0.5
                    xl = xlim;
                    xt = -0.5;
                    if isempty(find(xtks == xt))
                        idx_smaller = find(xtks < xt);
                        if ~isempty(idx_smaller) && xl(1) <= xt
                            xtks = [xtks(1:idx_smaller(end)) xt xtks(idx_smaller(end)+1:end)];
                            xtklabels = [xtklabels(1:idx_smaller(end)); {num2str(xt)}; xtklabels(idx_smaller(end)+1:end)];
                        end
                        clear idx_smaller
                    end
                    xt = 0.5;
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
                if xlim_t(1) < t_maxbef 
                    xlim_t(1) = t_maxbef;
                end
                if xlim_t(2) > t_maxaft 
                    xlim_t(2) = t_maxaft;
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
                        set(hl,'Position',[0.86 0.86 0.047 0.064],...
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
        
        
        trials_spectrogram_img = fullfile(savefolder, [animal '_' groupname '_trials_spect_' pdcond '_' dateofexp_yyyymmdd '_' bkstr '_trial' num2str(tri_str) '-' num2str(tri_end)]);
        saveas(gcf, trials_spectrogram_img, imFormat);
        
    end %plot_spectrogram

end
