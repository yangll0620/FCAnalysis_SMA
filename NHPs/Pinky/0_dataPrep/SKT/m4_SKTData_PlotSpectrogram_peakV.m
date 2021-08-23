function m4_SKTData_PlotSpectrogram_peakV()
%  extract lfp data respect to reachonset
% 
%  return:
%        lfptrials: nchns * ntemp * ntrials


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
animal = animal_extract(codecorresfolder);

%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');


savefig_format = 'tif';

% tdur used for calculate spectrogram respect to t_peakV
tdur_spect = [-0.6 1];


twin = 0.2;
toverlap = 0.18;
f_AOI = [8 40];

if strcmpi(animal, 'Kitty')
    t_AOI = [-0.3 0.3]; % tdur used for shown respect to t_peakV
end
if strcmpi(animal, 'Jo')
    t_AOI = [-0.3 0.2]; % tdur used for shown respect to t_peakV
end

%% save setup
savefolder = codecorresfolder;
savefilename = [animal 'peakV_timeStaComp'];

%% code start here

coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);



cond_cell = cond_cell_extract(animal);
for ci = 1: length(cond_cell)
    pdcond = cond_cell{ci};
    files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
    
    
    %%% extract psd_allchns_avgTrials:nf * nt * nchns
    psd_allchns_alltrials = []; 
    for fi = 1 : length(files)
        file = fullfile(files(fi).folder, files(fi).name);
        load(file, 'fs_lfp', 'lfpdata', 'goodTrials',  'T_idxevent_lfp',...
            'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma');
        
        
        nwin = round(twin * fs_lfp);
        noverlap = round(toverlap * fs_lfp);
        
        [~, ~, ntrials] = size(lfpdata);
        for tri = 1: ntrials
            
            % ignore trials marked with 0
            if ~goodTrials(tri)
                continue
            end
            
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            
            
            [~, idx] = max(smoothWspeed_trial(idx_reachonset_ma: idx_reach_ma, tri));
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
            
            % lfp_spect: lfp data used for calculating spectrogram (nchns * ntemp)
            lfp_spect = lfpdata(:, idx_spect_str:idx_spect_end, tri);
            
            % extract psd_allchns_plot using spectrogram, convert into dB and gaussfilt
            psd_allchns = []; % psd_allchns: nf * nt * nchns
            for chi = 1 : size(lfp_spect, 1)
                x = lfp_spect(chi, :);
                [~, freqs, times, psd] = spectrogram(x, nwin, noverlap,[],fs_lfp); % ps: nf * nt
                
                % convert into dB
                psd = 10 * log10(psd);
                
                % select freqs and corresponding psd
                idx_f = (freqs >= f_AOI(1) &  freqs <=f_AOI(2));
                freqs_plot =  freqs(idx_f);
                times = times + tdur_spect(1);
                idx_t = (times >= t_AOI(1) &  times <=t_AOI(2));
                times_plot_lfp = times(idx_t);
                
                psd_plot = psd(idx_f, :);
                psd_plot = psd_plot(:, idx_t);
                
                % cat into psd_allchns
                psd_allchns = cat(3, psd_allchns, psd_plot);
                
                clear x freqs times psd idx_f psd_plot idx_t
            end
            psd_allchns_alltrials = cat(4, psd_allchns_alltrials, psd_allchns);
            clear psd_allchns
            
            
            clear idx_reachonset_ma idx_reach_ma idx_peakV
            clear t_reachonset2peakV t_peakV2reach
        end
                
        clear file
        clear('fs_lfp', 'lfpdata', 'goodTrials', 'T_idxevent_lfp',...
            'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma');
        clear nchns ntrials tri
    end
    
    
    
    
    %%% plot psd_allchns_avgTrials:nf * nt * nchns
    load(fullfile(files(1).folder, files(1).name), 'T_chnsarea');
    
    % Group chns into STN, GP and others
    [clim_Spectrogram_STN, clim_Spectrogram_GP, clim_Spectrogram_Others] = clim_SKTSpectrogram_extract(animal);

    
    psd_allchns_avgTrials = mean(psd_allchns_alltrials, 4);
    for chi = 1 : size(psd_allchns_avgTrials, 3)
        brainarea = T_chnsarea.brainarea{chi};
        
        if contains(brainarea, 'stn')
            clim = clim_Spectrogram_STN;
        else
            if contains(brainarea, 'gp')
                clim = clim_Spectrogram_GP;
            else
                clim = clim_Spectrogram_Others;
            end
        end
        
        fig_sep = figure(); 
        set(fig_sep, 'PaperUnits', 'points',  'Position', [680 558 500 300]);
        ax_sep  = axes('Parent',fig_sep);
        
        psd_plot = psd_allchns_avgTrials(:, :,chi);
        psd_plot = imgaussfilt(psd_plot, 'FilterSize', [1 5]);
        imagesc(ax_sep, times_plot_lfp, freqs_plot, psd_plot);
        hold on
        
        colormap(jet)
        colorbar('FontSize', 9)
        
        ylabel('Frequency(Hz)', 'FontSize', 12, 'FontWeight', 'bold')
        xlabel('time/s', 'FontSize', 12, 'FontWeight', 'bold')
        xtls = xticklabels(ax_sep);
        xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char('peakV');
        xticklabels(ax_sep, xtls)
        set(ax_sep,'fontsize',11)
        set(ax_sep, 'Position', [0.09 0.15 0.8 0.79])
        
        title(ax_sep, [animal ' ' pdcond ':' brainarea])
        
        if isempty(clim)
            set(ax_sep,'YDir','normal')
        else
            set(ax_sep,'YDir','normal', 'CLim', clim)
        end
        
        % plot peakV line
        plot(ax_sep, [0 0], ylim, 'r--', 'LineWidth',1.5)
        
        % save
        savefile_sep = fullfile(savefolder, [animal '_peakV_' pdcond '_' brainarea]);
        saveas(fig_sep, savefile_sep, savefig_format);
        close(fig_sep)
            
        clear brainarea fig_sep ax_sep psd_plot
        clear xtls clim savefile_sep 
    end
    
    clear pdcond files
end











