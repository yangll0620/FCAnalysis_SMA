function m2_SKTData_SpectralAnalysis()
%% extract the Event Related Spectral perturbation (ERSP)



%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% input setup
animal = 'Pinky';

% Input folder
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_preprocessing');

% areas to be visulized
areas_vis = {'M1', 'lSMA', 'lVA', 'lVPLo','lCd', 'lVLo', 'rSMA', 'rMC','rVLo','rVPLo','rVA'};


t_trial = [-0.2 0.75];

% fft args
t_win = 0.25; t_overlap = t_win * 0.8;

%% save setup
savefolder = codecorresfolder;
savefilename_prefix = [animal '_STKData_ERSP'];


%% Starting Here
idxevent = table2array(T_idxevent);

n_win = round(t_win * fs); n_overlap = round(t_overlap * fs); nfft = n_win;

% calc and store ERSP in each area
for areai = 1: length(areas_vis)
    area = areas_vis{areai};
    
    % area chns
    chns = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
    
    % lfp1area: ntemporal * ntrials
    lfp1area = squeeze(mean(lfpdata(chns, :,:),1));
    
    for triali = 1: size(T_idxevent,2)
        idx_target = T_idxevent{triali,'TargetTime'};
        idx_mouth = T_idxevent{triali, 'MouthTimeix'};
        lfp1trial = lfp1area(idx_target:idx_mouth, triali);
    end
    
    % mean lfpdata(nchns * ntemporal * ntrials) across an area, lfp: 1 * ntemporal
    lfp = squeeze(mean(mean(lfpdata(chns, 500:end-250, :),1),3));
    
    % fft section
    yfft = fft(lfp);
    nfft = length(lfp);
    fs_fft = (0:nfft-1)*(fs/nfft);
    power = abs(yfft).^2/nfft;
    idx = find(f<=freqs(2) & f>=freqs(1));
    
    f_show = f(idx);
    power_show = power(idx);
    power_show = smooth(power_show, 20);
    
    plot(f_show, power_show);
    title(area);
    clear area chns lfp
end



