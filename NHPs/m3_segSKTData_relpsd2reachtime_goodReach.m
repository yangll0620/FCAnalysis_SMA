function m3_segSKTData_relpsd2reachtime_goodReach(animal, varargin)
%  spectrogram of all good Reach
%
%
%   Inputs:
%       animal
%
%       Name-Value: 
%           'F_AOI': frequency of Interested, default [8 40]
%
%           'codesavefolder' - code saved folder
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


% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
SKTSubfolder = 'SKT';
if strcmpi(animal, 'Kitty')
    SKTSubfolder = 'SKT_SegV';
end
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);



% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'f_AOI', [8 40], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
f_AOI = p.Results.f_AOI;

% copy code to savefolder if not empty
if ~isempty(codesavefolder)  
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

%%  input setup
image_type = 'tif';
% input folder
if strcmpi(animal, 'Kitty')
    SKTTrialfolder = 'm2_segSKTData_SelectTrials_goodReach';
end
if strcmpi(animal, 'Jo')
    SKTTrialfolder = 'm2_SKTData_SelectTrials';
end 
inputfolder = fullfile(codecorresParentfolder, SKTTrialfolder);


% align to event
align2 = SKTEvent.ReachOnset;
coli_align2 = uint32(align2);

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% starting: narrow filter the lfp data of all the files
cond_cell = cond_cell_extract(animal);
noisy_chns = noisy_chns_extract(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', codesavefolder);



% for each cond, extract reachtime and relative psd
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    if isempty(files)
        continue
    end
    
    eval(['tdur_trial = tdur_trial_' pdcond ';']);

    
    reachtime_1cond = [];
    relpsd_1cond = [];
    for fi = 1 : length(files)
        file = fullfile(files(fi).folder, files(fi).name);
        [reachtime_alltrials, relpsd_allchns_alltrials] = relpsd_eachTrials(file, tdur_trial, align2, noisy_chns);
        reachtime_1cond = cat(1, reachtime_1cond, reachtime_alltrials); 
        relpsd_1cond = cat(1, relpsd_1cond, relpsd_allchns_alltrials); 
        clear reachtime_alltrials relpsd_allchns_alltrials
        
        if ~exist('T_chnsarea', 'var')
            load(file, 'T_chnsarea');
            chnsOfI = chnOfInterest_extract(animal, 'codesavefolder', savecodefolder);
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);
        end
        
    end
    eval(['reachtime_' pdcond ' = reachtime_1cond;']);
    relpsd_1cond = relpsd_1cond(:, mask_chnOfI);
    eval(['relpsd_' pdcond ' = relpsd_1cond;']);
    
    clear pdcond files tdur_trial
    clear reachtime_1cond relpsd_1cond 
end

plotstyles = {'bo', 'k.'};
for chi = 1 : height(T_chnsarea)
    figure
    for i = 1 : length(cond_cell)
        pdcond = cond_cell{i};
        eval(['reachtimes = reachtime_' pdcond ';']);
        eval(['relpsds = relpsd_' pdcond '(:, chi);']);
        
        plot(relpsds, reachtimes, plotstyles{i}, 'DisplayName',pdcond);
        hold on
        
        clear reachtimes relpsds pdcond
    end
    title([animal ' ' T_chnsarea.brainarea{chi} ' psd vs reach time aligned to reachonset'])
    xlabel('relative psd')
    ylabel('reach time /s')
    xlim([-1 0.5])
    legend(cond_cell)
    savefile = [animal '_relpsd2reachtime_align2'  'reachonset_' T_chnsarea.brainarea{chi}];
    saveas(gcf, fullfile(savefolder, savefile), image_type);
end





% for each cond, extract reachtime and relative psd
tdur_spect = [-0.6 1]; % tdur used for calculate spectrogram respect to t_peakV
for i = 1 : length(cond_cell)
    pdcond = cond_cell{i};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    if isempty(files)
        continue
    end
    
    eval(['tdur_trial = tdur_trial_' pdcond ';']);

    
    reachtime_1cond = [];
    relpsd_1cond = [];
    for fi = 1 : length(files)
        file = fullfile(files(fi).folder, files(fi).name);
        [reachtime_alltrials, relpsd_allchns_alltrials] = relpsd_peakV_eachTrials(file, tdur_spect, noisy_chns);
        reachtime_1cond = cat(1, reachtime_1cond, reachtime_alltrials); 
        relpsd_1cond = cat(1, relpsd_1cond, relpsd_allchns_alltrials); 
        clear reachtime_alltrials relpsd_allchns_alltrials
        
        if ~exist('T_chnsarea', 'var')
            load(file, 'T_chnsarea');
            chnsOfI = chnOfInterest_extract(animal, 'codesavefolder', savecodefolder);
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);
        end
        
    end
    eval(['reachtime_' pdcond ' = reachtime_1cond;']);
    relpsd_1cond = relpsd_1cond(:, mask_chnOfI);
    eval(['relpsd_' pdcond ' = relpsd_1cond;']);
    
    clear pdcond files tdur_trial
    clear reachtime_1cond relpsd_1cond 
end

plotstyles = {'bo', 'k.'};
for chi = 1 : height(T_chnsarea)
    figure
    for i = 1 : length(cond_cell)
        pdcond = cond_cell{i};
        eval(['reachtimes = reachtime_' pdcond ';']);
        eval(['relpsds = relpsd_' pdcond '(:, chi);']);
        
        plot(relpsds, reachtimes, plotstyles{i}, 'DisplayName',pdcond);
        hold on
        
        clear reachtimes relpsds pdcond
    end
    title([animal ' ' T_chnsarea.brainarea{chi} ' psd vs reach time aligned to peakV'])
    xlabel('relative psd')
    ylabel('reach time /s')
    xlim([-1 0.5])
    legend(cond_cell)
    savefile = [animal '_relpsd2reachtime_align2'  'peakV_' T_chnsarea.brainarea{chi}];
    saveas(gcf, fullfile(savefolder, savefile), image_type);
end
close all


end



function [phtime_alltrials, relpsd_allchns_alltrials] = relpsd_eachTrials(file, tdur_trial, align2, noisy_chns, varargin)
% 
%   file: abs path for file
%   extract phase time and relative psd for each trial in file
%
%   Return:
%       phtime_alltrials: ntrials * 1
%
%       relpsd_allchns_alltrials: ntrials * nchns


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [15 20], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_AOI', [0 0.3], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_base', [-0.3 0], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 'phasetime', 'reachtime', @isstr);
parse(p,varargin{:});


f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
t_base = p.Results.t_base;
phasetime =  p.Results.phasetime;


coli_align2 = uint32(align2);


twin = 0.2;
toverlap = 0.15;

if strcmpi(phasetime, 'reachtime')
    i_strend = [2 3];
end



%% code start here
listOfVariables = who('-file', file);
if ismember('selectedTrials', listOfVariables)
    load(file,  'selectedTrials');
else
    if ismember('goodTrials', listOfVariables)
        load(file,  'goodTrials');
        selectedTrials = goodTrials;
        clear goodTrials
    end
end
clear listOfVariables


load(file, 'lfpdata', 'T_idxevent_lfp', 'fs_lfp','T_chnsarea');


% remove the noisy chns
mask_noisyChns = cellfun(@(x) contains(x, noisy_chns), T_chnsarea.brainarea);
mask_notDBS_notM1 = ~strcmp(T_chnsarea.brainarea, 'M1') & ~strcmp(T_chnsarea.electype, 'DBS');
mask_usedChns = ~(mask_noisyChns | mask_notDBS_notM1);
T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';

        
nwin = round(twin * fs_lfp);
noverlap = round(toverlap * fs_lfp);  


% extract reach time and relative psd for each trials
ntrials = length(selectedTrials);
phtime_alltrials = [];
relpsd_allchns_alltrials = [];
for tri = 1: ntrials
    
    % ignore trials marked with 0
    if ~selectedTrials(tri)
        continue
    end
    
    phtime = (T_idxevent_lfp{tri, i_strend(2)} -  T_idxevent_lfp{tri, i_strend(1)})/fs_lfp;
    
    %%% --- extract psd of lfp data for each chn: psd_allchns %%%    
    % extract trial with t_dur for lfp data: lfp_phase_1trial (nchns * ntemp)
    idxdur_lfp = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
    
    if iscell(lfpdata) % lfpdata: cell, lfpdata{tri}: nchns * ntemp
        tmp = lfpdata{tri}; % tmp: nchns * ntemp
    else 
        tmp = squeeze(lfpdata(:, :, tri)); % lfpdata: nchns * ntemp * ntrials
    end
    lfp_phase_1trial = tmp(mask_usedChns, idxdur_lfp(1):idxdur_lfp(2));
    clear idxdur_lfp  tmp
    
    % extract relative psd relpsd_allchns respect to t_base
    relpsd_allchns = []; % psd_allchns: nf * nt * nchns
    for chi = 1 : size(lfp_phase_1trial, 1)
        x = lfp_phase_1trial(chi, :);
        [~, freqs, times, psd] = spectrogram(x, nwin, noverlap,[],fs_lfp); % ps: nf * nt
        
        % convert into dB
        psd = 10 * log10(psd);
        
        % select freqs and corresponding psd
        idx_f_AOI = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        times = times + tdur_trial(1);
        idx_t_AOI = (times >= t_AOI(1) &  times <=t_AOI(2));
        idx_t_base = (times >= t_base(1) &  times <=t_base(2));
        
        psd_AOI = mean(mean(psd(idx_f_AOI, idx_t_AOI), 2),1);
        psd_base = mean(mean(psd(idx_f_AOI, idx_t_base), 2),1);
        relpsd = (psd_AOI - psd_base)/abs(psd_base);
        
        % cat into psd_allchns
        relpsd_allchns = cat(2, relpsd_allchns, relpsd);
        
        clear x freqs times psd 
        clear idx_f_AOI idx_t_AOI idx_t_base
        clear psd_AOI psd_base relpsd
    end

    phtime_alltrials = cat(1, phtime_alltrials, phtime);
    relpsd_allchns_alltrials = cat(1, relpsd_allchns_alltrials, relpsd_allchns);
    
    %%% final clear
    clear reachtime lfp_phase_1trial relpsd_allchns
end

end



function [phasetime_alltrials, relpsd_allchns_alltrials] = relpsd_peakV_eachTrials(file, tdur_spect, noisy_chns, varargin)
% 
%   file: abs path for file
%   extract phase time and relative psd for each trial in file
%
%   Return:
%       phasetime_alltrials: ntrials * 1
%
%       relpsd_allchns_alltrials: ntrials * nchns


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [15 20], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_AOI', [0 0.3], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_base', [-0.3 0], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 'phasetime', 'reachtime', @isstr);
parse(p,varargin{:});


f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
t_base = p.Results.t_base;
phasetime = p.Results.phasetime;


twin = 0.2;
toverlap = 0.15;

if strcmpi(phasetime, 'reachtime')
    i_strend = [2 3];
end



%% code start here
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);

load(file, 'lfpdata', 'T_idxevent_lfp', 'fs_lfp', 'selectedTrials', 'T_chnsarea',...
    'smoothWspeed_trial', 'Wrist_smooth_trial', 'T_idxevent_ma', 'fs_ma');


% remove the noisy chns
mask_noisyChns = cellfun(@(x) contains(x, noisy_chns), T_chnsarea.brainarea);
mask_notDBS_notM1 = ~strcmp(T_chnsarea.brainarea, 'M1') & ~strcmp(T_chnsarea.electype, 'DBS');
mask_usedChns = ~(mask_noisyChns | mask_notDBS_notM1);
T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';

        
nwin = round(twin * fs_lfp);
noverlap = round(toverlap * fs_lfp);  


% extract reach time and relative psd for each trials
ntrials = length(selectedTrials);
phasetime_alltrials = [];
relpsd_allchns_alltrials = [];
for tri = 1: ntrials
    
    % ignore trials marked with 0
    if ~selectedTrials(tri)
        continue
    end
    
    idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
    idx_reach_ma = T_idxevent_ma{tri, coli_reach};
    [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
    idx_peakV_ma = idx + idx_reachonset_ma -1;
    clear idx
    
    t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
    t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
    if t_reachonset2peakV < abs(t_AOI(1)) || t_peakV2reach < t_AOI(2)
        continue;
    end
    
    idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
    idx_spect_str = idx_peakV_lfp + tdur_spect(1) * fs_lfp;
    idx_spect_end = idx_peakV_lfp + tdur_spect(2) * fs_lfp;
    
   
    %%% --- extract psd of lfp data for each chn: psd_allchns %%%
    
    % lfp_spect: lfp data used for calculating spectrogram (nchns * ntemp)
    lfp_1trial = lfpdata{tri};
    lfp_spect = lfp_1trial(:, idx_spect_str:idx_spect_end);
    clear lfp_1trial

    % extract relative psd relpsd_allchns respect to t_base
    relpsd_allchns = []; % psd_allchns: nf * nt * nchns
    for chi = 1 : size(lfp_spect, 1)
        x = lfp_spect(chi, :);
        [~, freqs, times, psd] = spectrogram(x, nwin, noverlap,[],fs_lfp); % ps: nf * nt
        times = times + tdur_spect(1);
        
        % convert into dB
        psd = 10 * log10(psd);
        
        % select freqs and times
        idx_f_AOI = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        idx_t_AOI = (times >= t_AOI(1) &  times <=t_AOI(2));
        idx_t_base = (times >= t_base(1) &  times <=t_base(2));
        
        psd_AOI = mean(mean(psd(idx_f_AOI, idx_t_AOI), 2),1);
        psd_base = mean(mean(psd(idx_f_AOI, idx_t_base), 2),1);
        relpsd = (psd_AOI - psd_base)/abs(psd_base);
        
        % cat into psd_allchns
        relpsd_allchns = cat(2, relpsd_allchns, relpsd);
        
        clear x freqs times psd 
        clear idx_f_AOI idx_t_AOI idx_t_base
        clear psd_AOI psd_base relpsd
    end
    relpsd_allchns_alltrials = cat(1, relpsd_allchns_alltrials, relpsd_allchns);

    phtime = (T_idxevent_lfp{tri, i_strend(2)} -  T_idxevent_lfp{tri, i_strend(1)})/fs_lfp;
    phasetime_alltrials = cat(1, phasetime_alltrials, phtime);
    
    
    %%% final clear
    clear reachtime lfp_phase_1trial relpsd_allchns
end

end