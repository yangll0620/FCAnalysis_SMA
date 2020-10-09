%% Start Here 
if strcmp(eventName, 'reach')
    eventIdx_phase = reachIdx;
end

if strcmp(eventName, 'return')
    eventIdx_phase = returnIdx;
end

%  extract all lfp phase data
files = dir(fullfile(inputfolder, ['*_' cond '_*.mat']));
if strcmp(brainarea, 'M1')
    [lfps_phase, lfps_base, fs] = lfpm1_allfile(files, eventIdx_phase, tmin_phase, tmax_phase, tphase_bef, tbase_bef);
end
if strcmp(brainarea, 'SMA')
    [lfps_phase, lfps_base, fs] = lfpSMA_allfile(files, eventIdx_phase, tmin_phase, tmax_phase, tphase_bef, tbase_bef);
end
% calculate spectrogram
[spec_phase_psd, spec_phase_ts, spec_phase_fs] = spectrogram_calc(lfps_phase, twin, toverlap, fs, tphase_bef);
[spec_base_psd,  spec_base_ts, spec_base_fs] = spectrogram_calc(lfps_base, twin, toverlap, fs, tbase_bef);

% spec_phase_psd: nfs * nts * ntrials
psd_phase = spec_phase_psd;
psd_base = spec_base_psd;

% normalized psd_base and psd_phase
psd_total = sum(psd_phase, 1);
nfs = size(psd_phase, 1);
psd_phase_normalized = psd_phase ./ repmat(psd_total, nfs, 1, 1);
clear psd_total nfs

psd_total = sum(psd_base, 1);
nfs = size(psd_total, 1);
psd_base_normalized = psd_base ./ repmat(psd_total, nfs, 1, 1);
clear psd_total nfs

% mean psd across trials
psd_avgTrials = nanmean(psd_phase_normalized, 3);


% relative psd 
relativePSD = true;
if relativePSD
    psd_rel = psd_relative_calc(psd_phase_normalized, psd_base_normalized);
    
    psdrel_phase = nanmean(psd_rel, 3);
    
    clear psd_rel
end




%% plot
psd_phase = psd_avgTrials;


fig = figure;
sub1 = subplot(1,2,1,'Parent',fig);
sub2 = subplot(1,2,2, 'Parent', fig);


idx_F = spec_phase_fs> f_roi(1) & spec_phase_fs< f_roi(2);
idx_t = spec_phase_ts > t_roi(1);

fs_roi = spec_phase_fs(idx_F);
ts_roi = spec_phase_ts(idx_t) * 1000;


% ------ subplot spectrogram
psd_roi = psd_phase(idx_F, idx_t);


imagesc(sub1, ts_roi, fs_roi, psd_roi);
set(sub1,'YDir','normal') 

hold(sub1, 'on');
% plot zero line
plot(sub1, [0 0], sub1.YLim, 'r--')

%sub1.XTickLabel{1} = [eventName 'onset'];
xlabel(sub1,'ms')
ylabel(sub1, 'frequency')
title(sub1, [cond ' ' eventName ' in ' brainarea], 'FontSize', 12);
sub1.Position =  [0.08 0.1 0.75 0.8];
colorbar(sub1);



% ------ subplot relative psd
psdrel_roi = psdrel_phase(idx_F, idx_t);

idx_aft = ts_roi > 0;
psdrel_avgTs = mean(psdrel_roi(:, idx_aft), 2);

plot(sub2, psdrel_avgTs, fs_roi)
ylim(sub1.YLim);
sub2.Position =  [0.88 0.1 0.1 0.8];
xlabel(sub2, 'db')


%% save
saveas(gcf, ['psd_' brainarea '_' cond '_' eventName], 'png');
