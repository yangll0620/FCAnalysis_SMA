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

%% starting: narrow filter the lfp data of all the files
cond_cell = cond_cell_extract(animal);
noisy_chns = noisy_chns_extract(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', codesavefolder);
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
    end
    eval(['reachtime_' pdcond ' = reachtime_1cond;']);
    eval(['relpsd_' pdcond ' = relpsd_1cond;']);
    
    clear pdcond files tdur_trial
    clear reachtime_1cond reachtime_1cond
end

end



function [reachtime_alltrials, relpsd_allchns_alltrials] = relpsd_eachTrials(file, tdur_trial, align2, noisy_chns)
% 
%   file: abs path for file
%   extract reachtime and relative psd for each trial in file
%
%   Return:
%       reachtime_alltrials: ntrials * 1
%
%       relpsd_allchns_alltrials: ntrials * nchns


coli_align2 = uint32(align2);


twin = 0.2;
toverlap = 0.15;
f_AOI = [15 20];
t_AOI = [0 0.5];
t_base = [-0.5 0];


%% code start here
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
reachtime_alltrials = [];
relpsd_allchns_alltrials = [];
for tri = 1: ntrials
    
    % ignore trials marked with 0
    if ~selectedTrials(tri)
        continue
    end
    
    reachtime = (T_idxevent_lfp.TouchTimeix(tri) - T_idxevent_lfp.ReachTimeix(tri))/fs_lfp;
    
    %%% --- extract psd of lfp data for each chn: psd_allchns %%%    
    % extract trial with t_dur for lfp data: lfp_phase_1trial (nchns * ntemp)
    idxdur_lfp = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
    tmp = lfpdata{tri};
    lfp_phase_1trial = squeeze(tmp(:, idxdur_lfp(1):idxdur_lfp(2)));
    lfp_phase_1trial = lfp_phase_1trial(mask_usedChns, :, :);
    clear idxdur_lfp tmp
    
    % extract relative psd relpsd_allchns respect to t_base
    relpsd_allchns = []; % psd_allchns: nf * nt * nchns
    for chi = 1 : size(lfp_phase_1trial, 1)
        x = lfp_phase_1trial(chi, :);
        [~, freqs, times, psd] = spectrogram(x, nwin, noverlap,[],fs_lfp); % ps: nf * nt
        
        
        % select freqs and corresponding psd
        idx_f_AOI = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
        times = times + tdur_trial(1);
        idx_t_AOI = (times >= t_AOI(1) &  times <=t_AOI(2));
        idx_t_base = (times >= t_base(1) &  times <=t_base(2));
        
        psd_AOI = mean(mean(psd(idx_f_AOI, idx_t_AOI), 2),1);
        psd_base = mean(mean(psd(idx_f_AOI, idx_t_base), 2),1);
        relpsd = (psd_AOI - psd_base)/psd_base;
        
        % cat into psd_allchns
        relpsd_allchns = cat(2, relpsd_allchns, relpsd);
        
        clear x freqs times psd 
        clear idx_f_AOI idx_t_AOI idx_t_base
        clear psd_AOI psd_base relpsd
    end

    reachtime_alltrials = cat(1, reachtime_alltrials, reachtime);
    relpsd_allchns_alltrials = cat(1, relpsd_allchns_alltrials, relpsd_allchns);
    
    %%% final clear
    clear reachtime lfp_phase_1trial relpsd_allchns
end

end