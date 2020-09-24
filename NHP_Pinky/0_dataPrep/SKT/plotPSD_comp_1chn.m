function plotPSD_comp_1chn(lfp_normal, lfp_mild, lfp_moderate, F_range_roi, fs, savefolder, brainarea, task)

%%% extract psd_normal %%%
lfp  = lfp_normal;

% averaged across channels and all trials
lfp = squeeze(mean(mean(lfp, 3),1));

%
lfp = zscore(lfp);

% psd using pwelch
[psd, f_psd] = pwelch(lfp, length(lfp), 0, length(lfp), fs);

% normalized psd
psd = (psd - min(psd) )/ (max(psd) - min(psd));

psd_normal = psd;
clear psd lfp 


%%% extract psd_mild %%%
lfp  = lfp_mild;

% averaged across channels and all trials
lfp = squeeze(mean(mean(lfp, 3),1));

% psd using pwelch
[psd, ~] = pwelch(lfp, length(lfp), 0, length(lfp), fs);

% normalized psd
psd = (psd - min(psd) )/ (max(psd) - min(psd));

psd_mild = psd;
clear lfp psd



%%% extract psd_moderate %%%
lfp  = lfp_moderate;


% averaged across channels and all trials
lfp = squeeze(mean(mean(lfp, 3),1));

% zscore lfp
lfp = zscore(lfp);

% psd using pwelch
[psd, ~] = pwelch(lfp, length(lfp), 0, length(lfp), fs);

% normalized psd
psd = (psd - min(psd) )/ (max(psd) - min(psd));

psd_moderate = psd;
clear lfp psd




%%% --- plot ---%%%

% index for roi frequency
idx_roi = find(f_psd >= F_range_roi(1) & f_psd <= F_range_roi(2));

% fs_roi, and psd_roi 
fs_roi = f_psd(idx_roi);
psd_normal_roi = psd_normal(idx_roi);
psd_mild_roi = psd_mild(idx_roi);
psd_moderate_roi = psd_moderate(idx_roi);


smooth_span = 5;
plot(fs_roi, smooth(psd_normal_roi, smooth_span), 'DisplayName', 'normal')
hold on
plot(fs_roi, smooth(psd_mild_roi,smooth_span), 'DisplayName', 'mild')
plot(fs_roi, smooth(psd_moderate_roi, smooth_span), 'DisplayName', 'moderate')

legend()

% title
title([task ' PSD in ' brainarea ])

% save
% save
savefile = fullfile(savefolder, ['psd_' brainarea '_' task]);
saveas(gcf, savefile, 'pnd');