function [relpsd_alltrials] = relpsd_allTrials(file, chpair, varargin)
%    
%
%   file: abs path for file
%   extract phase time and relative psd for each trial in file
%
%   Inputs:
%       file
%       chnpair: 'stn1-2'
%
%       Name-Value: 
%           'F_AOI': frequency of Interested, default [8 40]
%           't_AOI': time duration of Interested, default [0 0.3]
%           't_base': time base for relative psd, default [-0.3 0]
%           'align2event': align to event {'reachOnset', 'peakV'}, i.e. time 0, default reachTime
%           'tdur_trial': time duration 
%
%   Return:
%
%       relpsd_alltrials: ntrials * 1, containing nan for bad trial


% parse params
p = inputParser;
addParameter(p, 'f_AOI', [15 20], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_AOI', [0 0.3], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 't_base', [-0.3 0], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 'align2event', 'reachOnset', @isstr);
addParameter(p, 'tdur_trial', [-0.5 0.5], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
parse(p,varargin{:});



f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
t_base = p.Results.t_base;
align2event =  p.Results.align2event; % align2event: 'reachOnset', 'peakV'
tdur_trial =  p.Results.tdur_trial;


twin = 0.2;
toverlap = 0.15;


%% code start here
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);

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


load(file, 'lfpdata', 'T_idxevent_lfp', 'fs_lfp','T_chnsarea', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');


% selected used chn
mask_usedChn = cellfun(@(x) strcmp(x,chpair), T_chnsarea.brainarea);

T_chnsarea = T_chnsarea(mask_usedChn, :); 
T_chnsarea.chni = [1: height(T_chnsarea)]';

        
nwin = round(twin * fs_lfp);
noverlap = round(toverlap * fs_lfp);  


% extract wanted phase time and relative psd for each trials
ntrials = length(selectedTrials);
relpsd_alltrials = [];
for tri = 1: ntrials
    
    % ignore trials marked with 0
    if ~selectedTrials(tri)
        relpsd_alltrials = cat(1, relpsd_alltrials, nan);
        continue
    end
    
    %%% --- extract psd of lfp data for each chn: psd_allchns %%%
    switch align2event
        case 'reachOnset'
            coli_align2 = uint32(SKTEvent.ReachOnset);
            idxdur_lfp = round(tdur_trial * fs_lfp) + T_idxevent_lfp{tri, coli_align2};
            clear coli_align2
        case 'peakV'
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            if iscell(smoothWspeed_trial) % smoothWspeed_trial: cell, lfpdata{tri}: ntemp * 1
                tmp = smoothWspeed_trial{tri}; % tmp: ntemp * 1
            else
                tmp = smoothWspeed_trial(:, tri); % smoothWspeed_trial: ntemp * ntrials
            end
            [~, idx] = max(tmp(idx_reachonset_ma: idx_reach_ma, 1));
            idx_peakV_ma = idx + idx_reachonset_ma -1;
            
            t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
            if t_reachonset2peakV < abs(t_AOI(1)) || t_peakV2reach < t_AOI(2)
                clear idx_reachonset_ma idx_reach_ma idx idx_peakV_ma t_reachonset2peakV t_peakV2reach tmp
                continue;
            end
            
            idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
            idxdur_lfp = round(tdur_trial * fs_lfp) + idx_peakV_lfp;
            clear idx_reachonset_ma idx_reach_ma idx idx_peakV_ma t_reachonset2peakV t_peakV2reach tmp
            clear idx_peakV_lfp
    end
    
    % extract trial with t_dur for lfp data: lfp_phase_1trial (nchns * ntemp)
    if iscell(lfpdata) % lfpdata: cell, lfpdata{tri}: nchns * ntemp
        tmp = lfpdata{tri}; % tmp: nchns * ntemp
    else 
        tmp = squeeze(lfpdata(:, :, tri)); % lfpdata: nchns * ntemp * ntrials
    end
    lfp_phase_1trial = tmp(mask_usedChn, idxdur_lfp(1):idxdur_lfp(2));
    clear idxdur_lfp  tmp
    
    % extract relative psd relpsd_allchns respect to t_base
    [~, freqs, times, psd] = spectrogram(lfp_phase_1trial, nwin, noverlap,[],fs_lfp); % ps: nf * nt
    
    % convert into dB
    psd = 10 * log10(psd);
    
    % select freqs and corresponding psd
    idx_f_AOI = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    times = times + tdur_trial(1);
    idx_t_AOI = (times >= t_AOI(1) &  times <=t_AOI(2));
    idx_t_base = (times >= t_base(1) &  times <=t_base(2));
    
    psd_AOI = mean(mean(psd(idx_f_AOI, idx_t_AOI), 2),1);
    psd_base = mean(mean(psd(idx_f_AOI, idx_t_base), 2),1);
    relpsd = (psd_AOI - psd_base)/abs(psd_base); % relpsd: scalar

    
    relpsd_alltrials = cat(1, relpsd_alltrials, relpsd);
    
    
    %%% final clear
    clear freqs times psd
    clear idx_f_AOI idx_t_AOI idx_t_base
    clear psd_AOI psd_base
    clear lfp_phase_1trial relpsd
end
plot(relpsd_alltrials, '*')
end