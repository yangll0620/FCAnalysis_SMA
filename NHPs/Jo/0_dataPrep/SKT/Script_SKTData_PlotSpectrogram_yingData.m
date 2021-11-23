clear
load('H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\root2\Ying Yu\SKB_EventRelatedBeta\Jo_data_SingleKluverBoardTask_Array_DBSLFP.mat')

tdur_trial = [-1 1];

dbslfp_trials_normal = [];
dbslfp_trials_mild = [];
dbslfp_trials_moderate = [];
for fi = 1 : length(data)
    pdstate  = data(fi).PDstate;
    
    dbslfp = data(fi).DBSlfp;
    fs_lfp = data(fi).lfpfs; 
    eventTime = data(fi).Eventtime;
    goodind = [1:size(eventTime, 1)];
    for gi = 1 : length(data(fi).goodind)
        goodind = intersect(goodind, data(fi).goodind{gi});
    end
    
    for gi = 1: length(goodind)
        tri = goodind(gi); 
        reachOnsetTime = eventTime(tri, 2);
        reachTime = eventTime(tri, 3);
        
        if (reachTime - reachOnsetTime) < 0.5 % skip if reach time < 0.5
            clear reachOnsetTime reachTime
            continue;
        end
        
        idx_str = round((reachOnsetTime + tdur_trial(1)) * fs_lfp);
        idx_end = round((reachOnsetTime + tdur_trial(2)) * fs_lfp);
        dbslfp_trial = dbslfp(idx_str: idx_end, :); % dbslfp_trial:  ntemp * nchns
        
        if strcmpi(pdstate, 'normal')
            dbslfp_trials_normal = cat(3, dbslfp_trials_normal, dbslfp_trial);
        else
            if strcmpi(pdstate, 'mild')
                dbslfp_trials_mild = cat(3, dbslfp_trials_mild, dbslfp_trial);
            else
                if strcmpi(pdstate, 'moderate')
                    dbslfp_trials_moderate = cat(3, dbslfp_trials_moderate, dbslfp_trial);
                end
            end
        end
        
        clear reachOnsetTime reachTime idx_str idx_end dbslfp_trial
    end
    
    if ~exist('fs', 'var')
        fs = fs_lfp;
    else
        if fs ~= fs_lfp
            disp(['fs in fi = ' num2str(fi) ' not consistent'])
        end
    end
    
    clear pdstate dbslfp fs_lfp eventTime tri
end
clear fi


%%% calculate spectrogram section
twin = 0.2;
toverlap = 0.18;
f_AOI = [8 40];
t_AOI = [-0.5 0.5];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);

clim = [-40 -10];

% calculate psd for each chn across trials
lfp_phase_trials = dbslfp_trials_normal;

psd_allchns = [];
for chi = 1 : size(lfp_phase_trials, 2)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfp_phase_trials, 3)
        x = lfp_phase_trials(:, chi, tri);
        x = zscore(x);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end
psd_allchns_normal = psd_allchns;
clear psd_allchns
clear lfp_phase_trials


% calculate psd for each chn across trials
lfp_phase_trials = dbslfp_trials_mild;

psd_allchns = [];
for chi = 1 : size(lfp_phase_trials, 2)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfp_phase_trials, 3)
        x = lfp_phase_trials(:, chi, tri);
        x = zscore(x);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end
psd_allchns_mild = psd_allchns;
clear psd_allchns
clear lfp_phase_trials



% calculate psd for each chn across trials
lfp_phase_trials = dbslfp_trials_moderate;

psd_allchns = [];
for chi = 1 : size(lfp_phase_trials, 2)
    psds = []; %  psds: nf * nt * ntrials
    for tri = 1: size(lfp_phase_trials, 3)
        x = lfp_phase_trials(:, chi, tri);
        x = zscore(x);
        [~, freqs, times, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    % convert into dB
    psds = 10 * log10(psds);
    psd = mean(psds, 3); % psd: nf * nt
    
    % select freqs/times and corresponding psd
    idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
    freqs_plot =  freqs(idx_f);
    
    times = times + tdur_trial(1);
    idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
    times_plot = times(idx_t);
    
    psd_plot = psd(idx_f, idx_t);
    
    % gauss filted
    psd_plot = imgaussfilt(psd_plot,'FilterSize',[3,11]);
    
    psd_allchns = cat(3, psd_allchns, psd_plot); % psd_allchns: nf * nt * nchns
    clear psds psd psd_plot idx_t idx_f
end
psd_allchns_moderate = psd_allchns;

clear psd_allchns
clear lfp_phase_trials

%%% plot spectrogram

for chi = 1 : 14
    if chi < 7
        chnname = ['stn' num2str(chi-1) '-' num2str(chi)];
    else
        chnname = ['gp' num2str(chi-1-7) '-' num2str(chi-7)];
    end
    
    % plot normal
    pdcond = 'normal';

    eval(['psd_allchns = psd_allchns_' pdcond ';'])
    fig_sep = figure();
    set(fig_sep, 'PaperUnits', 'points',  'Position', [680 558 600 300]);
    ax_sep  = axes('Parent',fig_sep);

    imagesc(ax_sep, times_plot, freqs_plot, psd_allchns(:, :, chi));
    colormap(jet)
    colorbar('FontSize', 9)
    ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
    xlabel('time/s', 'FontSize', 112, 'FontWeight', 'bold')
    xtls = xticklabels(ax_sep);
    xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = 'reachOnset';
    xticklabels(ax_sep, xtls)
    set(ax_sep,'fontsize',11)
    hold on
    % plot reach onset line
    plot(ax_sep, [0 0], ylim, 'r--', 'LineWidth',1.5)


    if isempty(clim)
        set(ax_sep,'YDir','normal')
    else
        set(ax_sep,'YDir','normal', 'CLim', clim)
    end
    title([chnname ': ' pdcond ])
    clear psd_allchns pdcond


    % plot mild
    pdcond = 'mild';

    eval(['psd_allchns = psd_allchns_' pdcond ';'])
    fig_sep = figure();
    set(fig_sep, 'PaperUnits', 'points',  'Position', [680 558 600 300]);
    ax_sep  = axes('Parent',fig_sep);

    imagesc(ax_sep, times_plot, freqs_plot, psd_allchns(:, :, chi));
    colormap(jet)
    colorbar('FontSize', 9)
    ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
    xlabel('time/s', 'FontSize', 112, 'FontWeight', 'bold')
    xtls = xticklabels(ax_sep);
    xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = 'reachOnset';
    xticklabels(ax_sep, xtls)
    set(ax_sep,'fontsize',11)
    hold on
    % plot reach onset line
    plot(ax_sep, [0 0], ylim, 'r--', 'LineWidth',1.5)


    if isempty(clim)
        set(ax_sep,'YDir','normal')
    else
        set(ax_sep,'YDir','normal', 'CLim', clim)
    end
    title([chnname ': ' pdcond ])
    clear psd_allchns pdcond


    % plot moderate
    pdcond = 'moderate';

    eval(['psd_allchns = psd_allchns_' pdcond ';'])
    fig_sep = figure();
    set(fig_sep, 'PaperUnits', 'points',  'Position', [680 558 600 300]);
    ax_sep  = axes('Parent',fig_sep);

    imagesc(ax_sep, times_plot, freqs_plot, psd_allchns(:, :, chi));
    colormap(jet)
    colorbar('FontSize', 9)
    ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
    xlabel('time/s', 'FontSize', 112, 'FontWeight', 'bold')
    xtls = xticklabels(ax_sep);
    xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = 'reachOnset';
    xticklabels(ax_sep, xtls)
    set(ax_sep,'fontsize',11)
    hold on
    % plot reach onset line
    plot(ax_sep, [0 0], ylim, 'r--', 'LineWidth',1.5)


    if isempty(clim)
        set(ax_sep,'YDir','normal')
    else
        set(ax_sep,'YDir','normal', 'CLim', clim)
    end
    title([chnname ': ' pdcond ])
    clear psd_allchns pdcond
end








