function m3_segSKTData_minMaxRelpsd_goodReach(animal, varargin)
%  min and max relative psd
%
% Kitty_TrialsWMarkers_moderate_20150408_bktdt2, trial = 8, 0  0  1  0
% Kitty_TrialsWMarkers_moderate_20150410_bktdt1, trial = 2, 1  0  0  0
% Kitty_TrialsWMarkers_moderate_20150417_bktdt1, trial = 9, 0  1  0  0
% Kitty_TrialsWMarkers_moderate_20150505_bktdt1, trial = 6, 0  0  0  1
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
pamars = inputParser;
addParameter(pamars, 'codesavefolder', '', @isstr);
addParameter(pamars, 'f_AOI', [15 20], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(pamars, 't_AOI', [0 0.3], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(pamars, 't_base', [-0.3 0], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);

parse(pamars,varargin{:});
codesavefolder = pamars.Results.codesavefolder;
f_AOI =  pamars.Results.f_AOI;
t_AOI =  pamars.Results.t_AOI;
t_base =  pamars.Results.t_base;

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


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

savefilename_prefix = 'relpsd2time_';


%% starting: narrow filter the lfp data of all the files
align2events = {'reachOnset'};
phasetimenames = {'reachTime'};

cond_cell = cond_cell_extract(animal);
noisy_chns = noisy_chns_extract(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', codesavefolder);

plotstyles = {'bo', 'r+', 'k*'};


%%% calculate phase time and relative psd
savematfile = fullfile(savefolder, [savefilename_prefix '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz' ...
                                    '_tAOI' num2str(t_AOI(1)) '-' num2str(t_AOI(2)) 's' ...
                                    '_tbase' num2str(t_base(1)) '-' num2str(t_base(2)) 's' '.mat']);
if ~exist(savematfile, 'file')
    for algi = 1 : length(align2events)
        align2event = align2events{algi};
        for phi = 1 : length(phasetimenames)
            phasetimename = phasetimenames{phi};
            
            disp(['align2event = ' align2event ', phasetimename = ' phasetimename])
            
            %%% for each cond, extract phase time and relative psd %%%
            for i = 1 : length(cond_cell)
                pdcond = cond_cell{i};
                
                files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
                
                if isempty(files)
                    continue
                end
                
                eval(['tdur_trial = tdur_trial_' pdcond ';']);
                
                
                phtime_1cond = [];
                relpsd_1cond = [];
                for fi = 1 : length(files)
                    file = fullfile(files(fi).folder, files(fi).name);
                    [phtime_alltrials, relpsd_allchns_alltrials] = find_relpsd(file, 'F_AOI', f_AOI, 't_AOI', t_AOI, 't_base', t_base, 'tdur_trial', tdur_trial, ...
                        'align2event', align2event, 'phasetimename', phasetimename, 'noisy_chns', noisy_chns);
                    phtime_1cond = cat(1, phtime_1cond, phtime_alltrials);
                    relpsd_1cond = cat(1, relpsd_1cond, relpsd_allchns_alltrials);
                    clear phtime_alltrials relpsd_allchns_alltrials
                    
                    if ~exist('T_chnsarea', 'var')
                        load(file, 'T_chnsarea');
                        chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
                        mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                        T_chnsarea = T_chnsarea(mask_chnOfI, :);
                    end
                    
                end
                eval(['phtime_' phasetimename '_' align2event '_' pdcond ' = phtime_1cond;']);
                relpsd_1cond = relpsd_1cond(:, mask_chnOfI);
                eval(['relpsd_' phasetimename '_' align2event '_' pdcond ' = relpsd_1cond;']);
                
                clear pdcond files tdur_trial
                clear phtime_1cond relpsd_1cond
            end
        end
    end
    save(savematfile, 'phtime_*', 'relpsd_*', 'T_chnsarea');
    clear('phtime_*', 'relpsd_*', 'T_chnsarea');
end

%%% 
load(savematfile, 'phtime_*', 'relpsd_*', 'T_chnsarea');
stnmask = cellfun(@(x) strcmp(x,'stn1-2'), T_chnsarea.brainarea);
relpsd_IOA = relpsd_reachTime_reachOnset_moderate(:, stnmask);
[sortpsd, indpsd] = sort(relpsd_IOA);
psdIOA = [sortpsd(1) sortpsd(2) sortpsd(end-1) sortpsd(end-2)]


pdcond = 'moderate';
eval(['tdur_trial = tdur_trial_' pdcond ';']);
for algi = 1 : length(align2events)
    align2event = align2events{algi};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    for fi = 1 : length(files)
        file = fullfile(files(fi).folder, files(fi).name);
        [phtime_alltrials, relpsd_allchns_alltrials] = find_relpsd(file, psdIOA, 'F_AOI', f_AOI, 't_AOI', t_AOI, 't_base', t_base, 'tdur_trial', tdur_trial, ...
            'align2event', align2event, 'noisy_chns', noisy_chns);
    end
end

end


function [phtime_alltrials, relpsd_allchns_alltrials] = find_relpsd(file, psdIOA, varargin)
%    
%
%   file: abs path for file
%   extract phase time and relative psd for each trial in file
%
%   Inputs:
%       file
%
%       Name-Value: 
%           'F_AOI': frequency of Interested, default [8 40]
%           't_AOI': time duration of Interested, default [0 0.3]
%           't_base': time base for relative psd, default [-0.3 0]
%           'align2event': align to event {'reachOnset', 'peakV'}, i.e. time 0, default reachTime
%           'tdur_trial': time duration 
%           'phasetimename': phase time name of interest {'reachTime', 'reachonset2PeakvTime', 'peakv2TouchTime'}, default 'reachTime'
%           'noisy_chns': noisy_chns, default {}
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
addParameter(p, 'align2event', 'reachOnset', @isstr);
addParameter(p, 'tdur_trial', [-0.5 0.5], @(x)isnumeric(x)&&isvector(x)&&length(x)==2);
addParameter(p, 'phasetimename', 'reachTime', @isstr);
addParameter(p, 'noisy_chns', {}, @iscell);
parse(p,varargin{:});


f_AOI = p.Results.f_AOI;
t_AOI = p.Results.t_AOI;
t_base = p.Results.t_base;
align2event =  p.Results.align2event; % align2event: 'reachOnset', 'peakV'
tdur_trial =  p.Results.tdur_trial;
phasetimename =  p.Results.phasetimename; % phasetimename: 'reachTime', 'reachonset2peakvTime', 'peakv2returnTime'
noisy_chns = p.Results.noisy_chns;


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


load(file, 'lfpdata', 'T_idxevent_lfp', 'fs_lfp','T_chnsarea', ...
    'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');


% only remain stn
stnmask = cellfun(@(x) strcmp(x,'stn1-2'), T_chnsarea.brainarea);
mask_usedChns = stnmask;
T_chnsarea = T_chnsarea(mask_usedChns, :);
        
nwin = round(twin * fs_lfp);
noverlap = round(toverlap * fs_lfp);  


% extract wanted phase time and relative psd for each trials
ntrials = length(selectedTrials);
phtime_alltrials = [];
relpsd_allchns_alltrials = [];
for tri = 1: ntrials
    
    % ignore trials marked with 0
    if ~selectedTrials(tri)
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
    lfp_phase_1trial = tmp(mask_usedChns, idxdur_lfp(1):idxdur_lfp(2));
    clear idxdur_lfp  tmp
    
    % extract relative psd relpsd_allchns respect to t_base
    x = lfp_phase_1trial;
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
    relpsd_1trial = (psd_AOI - psd_base)/abs(psd_base);
    
    if any(psdIOA == relpsd_1trial)
        [~,filename,~] = fileparts(file);
        disp([filename ', trial = ' num2str(tri) ', ' num2str(psdIOA == relpsd_1trial)])
        clear filename
    end
        
    clear x freqs times psd
    clear idx_f_AOI idx_t_AOI idx_t_base
    
    
    %%% final clear
    clear lfp_phase_1trial relpsd_1trial
end

end




