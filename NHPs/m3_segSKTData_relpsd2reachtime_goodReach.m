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


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% starting: narrow filter the lfp data of all the files
F_AOI = [15 20];
t_AOI = [0 0.3];
t_base = [-0.3 0];
align2events = {'reachOnset', 'peakV'};
phasetimenames = {'reachTime', 'reachonset2PeakvTime', 'peakv2TouchTime'};

cond_cell = cond_cell_extract(animal);
noisy_chns = noisy_chns_extract(animal);
[tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, 'codesavefolder', codesavefolder);

plotstyles = {'bo', 'r+', 'k*'};

for algi = 1 : length(align2events)
    align2event = align2events{algi};
    for phi = 1 : length(phasetimenames)
        phasetimename = phasetimenames{phi};
        
        disp(['align2event = ' align2event ', phasetimename = ' phasetimename])
        
        % for each cond, extract phase time and relative psd
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
                [phtime_alltrials, relpsd_allchns_alltrials] = relpsd_eachTrials(file, 'F_AOI', F_AOI, 't_AOI', t_AOI, 't_base', t_base, 'tdur_trial', tdur_trial, ...
                    'align2event', align2event, 'phasetimename', phasetimename, 'noisy_chns', noisy_chns);
                phtime_1cond = cat(1, phtime_1cond, phtime_alltrials);
                relpsd_1cond = cat(1, relpsd_1cond, relpsd_allchns_alltrials);
                clear reachtime_alltrials relpsd_allchns_alltrials
                
                if ~exist('T_chnsarea', 'var')
                    load(file, 'T_chnsarea');
                    chnsOfI = chnOfInterest_extract(animal, 'codesavefolder', savecodefolder);
                    mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
                    T_chnsarea = T_chnsarea(mask_chnOfI, :);
                end
                
            end
            eval(['phtime_' pdcond ' = phtime_1cond;']);
            relpsd_1cond = relpsd_1cond(:, mask_chnOfI);
            eval(['relpsd_' pdcond ' = relpsd_1cond;']);
            
            clear pdcond files tdur_trial
            clear phtime_1cond relpsd_1cond
        end
        
        % plot
        for chi = 1 : height(T_chnsarea)
            figure
            for i = 1 : length(cond_cell)
                pdcond = cond_cell{i};
                eval(['phtime = phtime_' pdcond ';']);
                eval(['relpsds = relpsd_' pdcond '(:, chi);']);
                
                if strcmpi(pdcond, 'normal')
                    pltstyle = plotstyles{1};
                else
                    if strcmpi(pdcond, 'mild')
                        pltstyle = plotstyles{2};
                    else
                        if strcmpi(pdcond, 'moderate')
                            pltstyle = plotstyles{3};
                        end
                    end
                end
                
                plot(relpsds, phtime, pltstyle, 'DisplayName',pdcond);
                hold on
                
                clear phtime relpsds pdcond pltstyle
            end
            title([animal ' ' T_chnsarea.brainarea{chi} ' psd vs ' phasetimename ' aligned to ' align2event])
            xlabel('relative psd')
            ylabel('reach time /s')
            xlim([-1 0.5])
            legend(cond_cell)
            
            % save
            savefile = [animal '_relpsd2' phasetimename '_align2' align2event '_' T_chnsarea.brainarea{chi}];
            saveas(gcf, fullfile(savefolder, savefile), image_type);
            close gcf
            
            clear savefile
        end
        
        clear phasetimename T_chnsarea phtime_*
    end
    clear align2event
end

end


function [phtime_alltrials, relpsd_allchns_alltrials] = relpsd_eachTrials(file, varargin)
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


% remove the noisy chns
mask_noisyChns = cellfun(@(x) contains(x, noisy_chns), T_chnsarea.brainarea);
mask_notDBS_notM1 = ~strcmp(T_chnsarea.brainarea, 'M1') & ~strcmp(T_chnsarea.electype, 'DBS');
mask_usedChns = ~(mask_noisyChns | mask_notDBS_notM1);
T_chnsarea = T_chnsarea(mask_usedChns, :); T_chnsarea.chni = [1: height(T_chnsarea)]';

        
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
    relpsd_allchns_alltrials = cat(1, relpsd_allchns_alltrials, relpsd_allchns);
    
    
    %%% extract phase time: phtime %%%     
    switch phasetimename
        case 'reachTime'
            phtime = (T_idxevent_ma{tri, coli_reach} -  T_idxevent_ma{tri, coli_reachonset})/fs_ma;
            
        case 'reachonset2PeakvTime'
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            if iscell(smoothWspeed_trial) % smoothWspeed_trial: cell, lfpdata{tri}: ntemp * 1
                tmp = smoothWspeed_trial{tri}; % tmp: ntemp * 1
            else
                tmp = smoothWspeed_trial(:, tri); % smoothWspeed_trial: ntemp * ntrials
            end
            [~, idx_peakV] = max(tmp(idx_reachonset_ma: idx_reach_ma, 1));
            phtime = idx_peakV / fs_ma;
            
            clear idx_reachonset_ma idx_reach_ma tmp idx_peakV
            
        case 'peakv2TouchTime'
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            if iscell(smoothWspeed_trial) % smoothWspeed_trial: cell, lfpdata{tri}: ntemp * 1
                tmp = smoothWspeed_trial{tri}; % tmp: ntemp * 1
            else
                tmp = smoothWspeed_trial(:, tri); % smoothWspeed_trial: ntemp * ntrials
            end
            [~, idx_peakV] = max(tmp(idx_reachonset_ma: idx_reach_ma, 1));
            phtime = (length(tmp) - idx_peakV) / fs_ma;
            
            clear idx_reachonset_ma idx_reach_ma tmp idx_peakV    
    end
    phtime_alltrials = cat(1, phtime_alltrials, phtime);
    
    
    %%% final clear
    clear reachtime lfp_phase_1trial relpsd_allchns
end

end



