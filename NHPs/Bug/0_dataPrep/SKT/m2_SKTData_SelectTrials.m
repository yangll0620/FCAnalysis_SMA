function m2_SKTData_SelectTrials()
%  plot spectrogram of each trial and manually select trials
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
    
    tdur_trial_normal = [-0.6 2];
    tdur_trial_moderate = [-0.6 2];
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

% align to event
align2 = SKTEvent.ReachOnset;
coli_align2 = uint32(align2);

%% save setup
savefolder = codecorresfolder;


%% starting
files = dir(fullfile(inputfolder, '*.mat'));
filei = 19;
nfiles = length(files);
while(filei <=  nfiles)
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent', 'fs', 'T_chnsarea');
    
    
    % extract dateofexp, bktdt and pdcond
    idx_cond = cellfun(@(x) contains(filename, x), cond_cell);
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
    
    [~, ~, ntrials] = size(lfpdata);
    savefile_goodTrialsMarkers = fullfile(savefolder, [animal '_TrialsWMarkers_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_' bkstr  '.mat']);
    
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
    goodTrials_Loaded = false;
    if exist(savefile_goodTrialsMarkers, 'file')
        load(savefile_goodTrialsMarkers, 'tbl_goodTrialsMarks', 'idxGroups');
        groupNames_load = tbl_goodTrialsMarks.Properties.VariableNames;
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
            
            goodTrials_Loaded = true;
        end
    end
    idxGroups = idxGroups_split;
    groupNames = groupNames_split;
    if goodTrials_Loaded
        goodTrials_allGs = tbl_goodTrialsMarks{:, :};
    else
        goodTrials_allGs = ones(ntrials, length(idxGroups));
    end
    clear idxGroups_split groupNames_split
    
    showname = ['fi = ' num2str(filei) ':' animal '-' pdcond '-' datestr(dateofexp, 'yyyymmdd') '-' bkstr];
    
    % check_spectrogram_oneGroup
    [clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);
    for idxGi = 1 : length(idxGroups)
        idxs = idxGroups{idxGi};
        lfpdata_1group = lfpdata(idxs, :, :);
        T_chnsarea_1group = T_chnsarea(idxs, :);
        
        if all(contains(T_chnsarea_1group.brainarea, 'stn'))
            clim_Spectrogram = clim_Spectrogram_STN;
        else
            if all(contains(T_chnsarea_1group.brainarea, 'gp'))
                clim_Spectrogram = clim_Spectrogram_GP;
            else
                clim_Spectrogram = clim_Spectrogram_Others;
            end
        end
        
        goodTrials_1Grp = goodTrials_allGs(:, idxGi);
        goodTrials_1Grp = check_spectrogram_oneGroup(lfpdata_1group, T_idxevent, T_chnsarea_1group, fs, goodTrials_1Grp, showname, clim_Spectrogram);
        
        savefolder_eachTrial = fullfile(savefolder, 'eachTrial');
        if ~exist(savefolder_eachTrial, 'dir')
            mkdir(savefolder_eachTrial);
        end
        
        goodTrials_allGs(:, idxGi) = goodTrials_1Grp;
        
        clear idxs lfpdata_1group T_chnsarea_1group goodTrials_1Grp clim_Spectrogram
    end
    
    % goodTrials = goodTrials_allGs(:, 1) & goodTrials_allGs(:, 2) & goodTrials_allGs(:, 3)
    goodTrials = ones(size(goodTrials_allGs, 1), 1);
    for gi = 1 : size(goodTrials_allGs, 2)
        goodTrials = goodTrials & goodTrials_allGs(:, gi);
    end
    tbl_goodTrialsMarks = array2table(goodTrials_allGs, 'VariableNames', groupNames);
    save(savefile_goodTrialsMarkers, 'goodTrials', 'tbl_goodTrialsMarks', 'idxGroups', 'lfpdata', 'T_idxevent', 'fs', 'T_chnsarea');
    
    
    %%% --- show one day spectrogram using goodTrials --- %%%
    if strcmp(pdcond, 'normal')
        tdur_trial = tdur_trial_normal;
        t_minmax_reach = t_minmax_reach_normal;
        t_minmax_return = t_minmax_return_normal;
    end
    if strcmp(pdcond, 'mild')
        tdur_trial = tdur_trial_mild;
        t_minmax_reach = t_minmax_reach_mild;
        t_minmax_return = t_minmax_return_mild;
    end
    if strcmp(pdcond, 'moderate')
        tdur_trial = tdur_trial_moderate;
        t_minmax_reach = t_minmax_reach_moderate;
        t_minmax_return = t_minmax_return_moderate;
    end
    
    
    % extract phase trials respected to coli_align2
    lfp_phase_trials = [];
    for tri = 1: size(lfpdata, 3)
        
        % ignore trials marked with 0
        if ~goodTrials(tri)
            continue
        end
        
        % select trials based on reach and return duration
        t_reach = (T_idxevent{tri, coli_reach} - T_idxevent{tri, coli_reachonset}) / fs;
        t_return = (T_idxevent{tri, coli_mouth} - T_idxevent{tri, coli_returnonset}) / fs;
        if t_reach < t_minmax_reach(1) || t_reach > t_minmax_reach(2)
            goodTrials(tri) = 0;
            clear t_reach
            continue
        end
        if t_return < t_minmax_return(1) || t_reach > t_minmax_return(2)
            goodTrials(tri) = 0;
            clear t_return
            continue
        end
        
        % extract trial with t_dur
        idxdur = round(tdur_trial * fs) + T_idxevent{tri, coli_align2};
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri); % lfp_phase_1trial: nchns * ntemp
        
        lfp_phase_trials = cat(3, lfp_phase_trials, lfp_phase_1trial); % lfp_phase_trials: nchns * ntemp * ntrials
        
        clear t_reach t_return idxdur lfp_phase_1trial
    end
    
    if isempty(lfp_phase_trials)
        disp('lfp_phase_trials is empty, skip spectrogram across trials')
        filei = filei + 1;
        continue;
    end
    
    % plot spectrogram across all good trials of one day
    disp(['# of Good Trials for averaged spectrogram across trials: ' num2str(size(lfp_phase_trials, 3))])
    idxGroupNames = tbl_goodTrialsMarks.Properties.VariableNames;
    plot_spectrogram_acrossTrials(lfp_phase_trials, T_chnsarea, idxGroups, idxGroupNames, tdur_trial, fs, animal, pdcond, align2, showname)
    clear idxGroupNames
    
    % Recheck today or Check the next day
    reply = input(['Check the next day (y) or Recheck this day (n)[y]:'], 's');
    if isempty(reply) || lower(reply) == 'y'
        saveas(gcf, fullfile(savefolder, [animal '_goodTrials_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_' bkstr]), 'png');
        
        filei = filei + 1;
    end
    close all
    
    
    clear filename lfpdata T_idxevent fs T_chnsarea tbl_goodTrialsMarks
    clear mask_STN mask_GP mask_Others idxGroups
    clear idx_cond datestrmatch pdcond dateofexp bktdt
    clear tdur_trial t_minmax_reach t_minmax_return
    clear lfpdata_origin T_idxevent_origin showname
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
fig_left = 50;
fig_bottom = 50;
fig_width = 1800;
fig_height = 900;


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


function goodTrials = check_spectrogram_oneGroup(lfpdata, T_idxevent, T_chnsarea, fs, goodTrials, showname, clim)
% lfpdata: nchns * ntemp * ntrials


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

% subplot/Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 1800;
fig_height = 900;


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


% trial number checkbox parameters
cb_width = 150;
cb_height = 50;


[nchns, ~, ntrials] = size(lfpdata);
ngs = ceil(ntrials / subp_ntrials);

fig = figure(); set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
hcbs = zeros(ntrials, 1);
checkedAllGs = zeros(ngs, 1);


gi = 0;
clf(fig);
c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_btn_next_left pos_btn_next_bottom btn_width btn_height]);
c_next.Callback = @nextTrials;
title(['gi = ' num2str(gi + 1) '/' num2str(ngs)])
tri_str = 1;
tri_end = subp_ntrials;
plot_spectrogram()
checkedAllGs(gi + 1) = 1;

uiwait(fig)


    function showNewTrials()
        clf(fig);
        
        if(gi < ngs-1)
            c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_btn_next_left pos_btn_next_bottom btn_width btn_height]);
            c_next.Callback = @nextTrials;
        end
        
        if(gi > 0)
            c_prev = uicontrol(fig, 'Style','pushbutton', 'String', 'Previous', 'Position', [pos_btn_prev_left pos_btn_prev_bottom btn_width btn_height]);
            c_prev.Callback = @prevTrials;
        end
        
        % start and end trial number for next graph
        tri_str = gi * subp_ntrials + 1;
        tri_end = (gi + 1) * subp_ntrials;
        if tri_end > ntrials
            tri_end = ntrials;
        end
        plot_spectrogram()
        
        checkedAllGs(gi + 1) = 1;
        if all(checkedAllGs)
            c_Finish = uicontrol(fig, 'Style','pushbutton', 'String', 'Finish', 'Position', [pos_btn_finish_left pos_btn_finish_bottom btn_width btn_height]);
            c_Finish.Callback = @finishCheck;
        end
        
    end

    function finishCheck(src,event)
        finishedCheck = true;
        close(fig)
    end % finishCheck

    function nextTrials(src,event)
        gi = mod(gi + 1, ngs);
        showNewTrials()
    end % nextTrials

    function prevTrials(src,event)
        gi = mod(gi - 1, ngs);
        showNewTrials()
    end % prevTrials

    function box_value(hObj,event)
        % Called when boxes are used
        
        % find the used box
        idx = find(hcbs==hObj);
        
        goodTrials(idx) =  get(hObj,'Value');
    end

    function plot_spectrogram()
        % plot lfpdata_1group of all the channels: nchns * ntemp * ntrial
        
        annotation(gcf,'textbox',...
            [subp_startLeft subp_startTop 1 0.03],...
            'String', {showname}, ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
        
        for tri = tri_str: tri_end
            coli = mod(tri,subp_ntrials);
            if coli == 0
                coli = subp_ntrials;
            end
            
            % left postion for tri
            subp_left = (coli -1) * (subp_width + supb_deltaX )+ subp_startLeft;
            
            % plot trial number Name and Checkbox
            pos_cb_trialN_left = (subp_left + subp_width/2) * fig_width;
            pos_cb_trialN_bottom = (subp_startTop) * fig_height;
            hcbs(tri) = uicontrol('Style','checkbox','Value', goodTrials(tri),...
                'Position', [pos_cb_trialN_left pos_cb_trialN_bottom cb_width cb_height],'String', ['triali = ' num2str(tri) '/' num2str(ntrials)]);
            set(hcbs(tri),'Callback',{@box_value});
            
            
            for chi = 1 : nchns
                
                x = lfpdata(chi, 1:T_idxevent{tri, coli_mouth}, tri);
                
                % spectrogram
                [~, freqs, times, psd] = spectrogram(x, round(twin * fs), round(toverlap * fs),[],fs); % psd: nf * nt
                
                % select freqs_plot and corresponding psd
                idx_f = (freqs >= f_AOI(1) & freqs <= f_AOI(2));
                freqs_plot =  freqs(idx_f);
                psd_plot = psd(idx_f, :);
                % convert to db and gaussfilt
                psd_plot = 10 * log10(psd_plot);
                psd_plot = imgaussfilt(psd_plot,'FilterSize',5);
                times_plot = times - T_idxevent{tri, coli_align2} / fs;
                
                
                % spectrogram subplot
                subp_bottom = subp_startTop - subp_height - (chi -1) * (subp_height + supb_deltaY);
                subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
                imagesc(times_plot, freqs_plot, psd_plot);
                set(gca,'YDir','normal', 'CLim', clim)
                colormap(jet)
                colorbar
                if (chi == nchns)
                    xlabel('time/s')
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
                for eventi = 1 : width(T_idxevent)
                    t_event = (T_idxevent{tri, eventi} - T_idxevent{tri, coli_align2}) / fs;
                    plot([t_event t_event], ylim, '--', 'LineWidth',1.5)
                    clear t_event
                end
                clear eventi
                
                
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
        end
        
        
    end %plot_spectrogram

end
