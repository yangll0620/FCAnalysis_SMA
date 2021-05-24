function m5_imCohUsingFFT_EventPhase()
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

% animal
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '/', '[A-Za-z]*']);
elseif ispc
    % Code to run on Windows platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '\\', '[A-Za-z]*']);
else
    disp('Platform not supported')
end
animal = codecorresfolder(fi + length('NHPs') + 1:j);

%% save setup
savefolder = codecorresfolder;


%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');

fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;

EventPhases = {'preMove'; 'Reach'; 'Return'};

image_type = 'bmp';


twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];

cond_cell = cond_cell_extract(animal);
[t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = ...
    goodSKTTrials_reachReturn_tcritiria(animal);
if strcmpi(animal, 'bug')
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
end
if strcmpi(animal, 'jo')
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
end

if strcmpi(animal, 'pinky')
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
end

for ei = 1: length(EventPhases)
    event = EventPhases{ei};
    if strcmpi(event, 'preMove')
        align2 = SKTEvent.TargetOnset;
        t_AOI = [-0.2 0];
    end
    if strcmpi(event, 'Reach')
        align2 = SKTEvent.ReachOnset;
        t_AOI = [0 0.2];
    end
    if strcmpi(event, 'Return')
        align2 = SKTEvent.ReturnOnset;
        t_AOI = [0 0.2];
    end
    
    
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        savefile_prefix = fullfile(savefolder, [animal '_' pdcond '_' event '_align2' char(align2)]);
        
        if ~exist([savefile_prefix, '.mat'])
            if strcmp(pdcond, 'normal')
                t_minmax_reach = t_minmax_reach_normal;
                t_minmax_return = t_minmax_return_normal;
                tdur_trial = tdur_trial_normal;
            else
                if strcmp(pdcond, 'mild')
                    t_minmax_reach = t_minmax_reach_mild;
                    t_minmax_return = t_minmax_return_mild;
                    tdur_trial = tdur_trial_mild;
                else
                    if strcmp(pdcond, 'moderate')
                        t_minmax_reach = t_minmax_reach_moderate;
                        t_minmax_return = t_minmax_return_moderate;
                        tdur_trial = tdur_trial_moderate;
                    end
                end
                
            end
            
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end
            
            [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach, t_minmax_return);
            
            [nchns, ~, ntrials] = size(lfptrials);
            
            % calculate imaginary of Coherency
            for chni = 1 : nchns-1
                lfptrialsi = squeeze(lfptrials(chni, :, :));
                for chnj = chni : nchns
                    lfptrialsj = squeeze(lfptrials(chnj, :, :));
                    
                    cross_density_sum = 0;
                    for triali = 1: ntrials
                        x = lfptrialsi(: , triali);
                        y = lfptrialsj(: , triali);
                        
                        [Sx, fx, tx, ~] = spectrogram(x, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
                        [Sy, ~, ~, ~] = spectrogram(y, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
                        
                        if chni == 1 && chnj == chni && triali == 1
                            freqs = fx;
                            times = tx + tdur_trial(1);
                            idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
                            idx_t = (times>=t_AOI(1) &  times<=t_AOI(2));
                            f_selected = round(freqs(idx_f), 3);
                            clear freqs times
                        end
                        
                        phix = angle(Sx(idx_f, idx_t));
                        phiy = angle(Sy(idx_f, idx_t));
                        cross_density_sum = cross_density_sum + exp(1i * (phix - phiy));
                        
                        clear x y Sx fx tx Sy
                        clear phix phiy
                        
                    end
                    
                    if chni == 1 && chnj == chni + 1
                        [nf, nt] = size(cross_density_sum);
                        cross_densitys = zeros(nchns, nchns, nf, nt);
                        clear nf nt
                    end
                    cross_densitys(chni, chnj, :, :) = cross_density_sum / ntrials; % cross_densitys: nchns * nchns * nf  * nt
                    cross_densitys(chnj, chni, :, :) = cross_densitys(chni, chnj, :, :);
                    
                    
                    clear cross_density_sum lfptrialsj
                end
                
                clear lfptrialsi
            end
            iCoh = imag(cross_densitys);
            iCoh_acrossTime = mean(iCoh, 4);
            
            lfp1 = squeeze(lfptrials(1, :, :));
            lfp2 = squeeze(lfptrials(10, :, :));
            [mus, stds] = psedo_SKTLFP_Test_acrossTime(lfp1, lfp2, 100, fs, twin, toverlap, f_AOI, t_AOI - tdur_trial(1));
            clear lfp1 lfp2
            
            
            % pvalues using permutation test
            [nchns, ~, nf] = size(iCoh_acrossTime);
            pvals = zeros(size(iCoh_acrossTime));
            for fi = 1 : nf
                
                mu = mus(fi,1);
                std = stds(fi, 1);
                pd = makedist('Normal','mu',mu,'sigma',std);
                
                for chni = 1: nchns -1
                    for chnj = chni : nchns
                        x = iCoh_acrossTime(chni, chnj, fi);
                        pvals(chni, chnj, fi) = (1 - cdf(pd,abs(x)));
                        pvals(chnj, chni, fi) = pvals(chni, chnj, fi);
                        clear x
                    end
                end
                
                clear mu std pd
            end
            % Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR)
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
            
            % set values not significant as 0
            iCoh_acrossTime(h==0) = 0;
            
            % show and save iCoh images
            
            % generate chnPairNames, such as M1-stn0-1
            nf = size(iCoh_acrossTime, 3);
            chnPairNames = {};
            iCoh_1time = zeros(nchns * (nchns -1)/2, nf);
            ci = 0;
            for chni = 1 : nchns -1
                for chnj = chni + 1  : nchns
                    chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
                    
                    ci = ci + 1;
                    iCoh_1time(ci, :) = iCoh_acrossTime(chni, chnj, :);
                end
            end
            
            % save data
            save([savefile_prefix '.mat'], 'iCoh_1time', 'f_selected', 'chnPairNames', 'ntrials')
            
        else
            load([savefile_prefix '.mat'], 'iCoh_1time', 'f_selected',  'chnPairNames', 'ntrials')
        end
        
        
        
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames); 
        usedChnPairsMask = M1DBS_mask | STN2GP_mask;
        clear M1DBS_mask STN2GP_mask
        
        showData = abs(iCoh_1time(usedChnPairsMask, :));
        chnPairNames_show = chnPairNames(usedChnPairsMask);
        
        % plot
        figure;
        set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
        imagesc(showData)
        set(gca, 'Position', [0.09 0.05 0.9 0.9])
        [npairs, nf] = size(showData);
        xticks([1:nf])
        xticklabels(round(f_selected,2))
        yticks([1:npairs]);
        set(gca,'YTickLabel',chnPairNames_show,'fontsize',12,'FontWeight','bold')
        xlabel('freqs')
        title([animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' 'align2 = ' char(align2) 'ntrials = ' num2str(ntrials)])
        set(gca,'CLim', [0 1])
        colorbar
        
        
        
        chnPair_prev = '';
        for ci = 1: length(chnPairNames_show)
            chnPair = chnPairNames_show{ci};
            
            % replace M1-stn0-1 to M1-STN
            s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
            if ~isempty(s_stn)
                for si = 1 : length(s_stn)
                    chnPair = strrep(chnPair, s_stn{si}, 'STN');
                end
            end
            % replace M1-stn0-1 to M1-STN
            s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
            if ~isempty(s_gp)
                for si = 1 : length(s_gp)
                    chnPair = strrep(chnPair, s_gp{si}, 'GP');
                end
            end
            
            if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
                hold on; plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'w--')
                % Create line
            end
            chnPair_prev = chnPair;
            
            clear s_stn s_gp chnPair
        end
        
        % save image
        saveas(gcf, fullfile(savefolder, [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type]), image_type);
        close all
        
        
        clear pdcond t_minmax_reach t_minmax_return tdur_trial
        clear files lfptrials fs T_chnsarea nchns ntrials
        clear iCoh iCoh_acrossTime mus stds
        clear idx_f idx_t f_selected t_selected
        clear nf chnPairNames iCoh_1time ci
        clear M1DBS_mask STN2GP_mask usedChnPairsMask showData npairs
        clear savefile_prefix
    end
    
    
    % combined all pdconds
    images_combined = [];
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        try
            image = imread(fullfile(savefolder, [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type]));
            [m, n, c] =  size(image);
        catch
            if exist('m', 'var')
                image = zeros(m,n,c) + 255;
            end
        end
        
        images_combined = cat(2, images_combined, image);
        imwrite(images_combined,fullfile(savefolder, [animal '_conn_combined_' event '_align2' char(align2) '.' image_type]))
    end
    
    
    clear t_AOI align2 event
end




%%%-------  plot Peak ---------%
patterns = {'ko', 'g+', 'r*'};
cond_cell = cond_cell_extract(animal);

for ei = 1: length(EventPhases)
    event = EventPhases{ei};
    if strcmpi(event, 'preMove')
        align2 = SKTEvent.TargetOnset;
    end
    if strcmpi(event, 'Reach')
        align2 = SKTEvent.ReachOnset;
    end
    if strcmpi(event, 'Return')
        align2 = SKTEvent.ReturnOnset;
    end
    
    figure
    set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
    for ci = 1 : length(cond_cell)
        pdcond = cond_cell{ci};
        savefile = fullfile(savefolder, [animal '_' pdcond '_' event '_align2' char(align2) '.mat']);
        
        load(savefile, 'chnPairNames', 'iCoh_1time');
        if ~exist('f_selected_show', 'var')
            load(savefile, 'f_selected');
            f_selected_show = f_selected;
        else
            load(savefile, 'f_selected');
            if any(~(f_selected_show == f_selected))
                disp([pdcond ' f_selected_show not equal f_selected'])
                continue;
            end
            clear f_selected
        end



        % select used chns
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        usedChnPairsMask = M1DBS_mask | STN2GP_mask;
        showData = abs(iCoh_1time(usedChnPairsMask, :));
        if ~exist('chnPairNames_show', 'var')
            chnPairNames_show = chnPairNames(usedChnPairsMask);
        else
            chnPairNames_show_new = chnPairNames(usedChnPairsMask);
            cellcmp = strcmp(chnPairNames_show, chnPairNames_show_new);
            if any(~cellcmp)
                disp([pdcond ' chnPairNames_show_new not equal chnPairNames_show'])
                continue;
            end
            clear chnPairNames_show_new
        end

        if ~exist('npairs', 'var')
            [npairs, nf] = size(showData);
        else
            [m, n] = size(showData);
            if m ~= npairs || n ~= nf 
                disp([pdcond ' m ~= npairs or n ~= nf'])
                continue;
            end
            clear m n
        end

        [iCoh_Peak, idxs] = max(abs(showData), [], 2);
        idxs = idxs(find(iCoh_Peak ~= 0));
        pairs = find(iCoh_Peak ~= 0);

        % plot
        plot(f_selected_show(idxs), pairs, patterns{ci}, 'DisplayName',pdcond)
        ylim([1, npairs])
        xlim([f_selected_show(1), f_selected_show(end)])
        hold on



        if ci == length(cond_cell)
            chnPair_prev = '';
            for ci = 1: length(chnPairNames_show)
                chnPair = chnPairNames_show{ci};

                % replace M1-stn0-1 to M1-STN
                s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
                if ~isempty(s_stn)
                    for si = 1 : length(s_stn)
                        chnPair = strrep(chnPair, s_stn{si}, 'STN');
                    end
                end
                % replace M1-stn0-1 to M1-STN
                s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
                if ~isempty(s_gp)
                    for si = 1 : length(s_gp)
                        chnPair = strrep(chnPair, s_gp{si}, 'GP');
                    end
                end

                if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
                    hold on; 
                    plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'k--')
                    % Create line
                end
                chnPair_prev = chnPair;

                clear s_stn s_gp chnPair
            end
        end


        clear M1DBS_mask STN2GP_mask usedChnPairsMask
        clear chnPairNames f_selected iCoh_1time
    end
    yticks([1:npairs])
    set(gca,'YTickLabel',chnPairNames_show,'fontsize',12,'FontWeight','bold')
    set(gca, 'YDir','reverse')
    xticks(f_selected_show)
    xticklabels(round(f_selected_show,2))
    set(gca, 'Color', 'w');
    hleg = legend('show', 'Color', 'w');
    hleg.String(end-1:end) = [];
    set(gca, 'Position', [0.09 0.05 0.9 0.9])
    title([ animal ' Peak connectivity in ' event])

    % save png
    saveas(gcf, fullfile(savefolder,[animal '_peakFC_' event '.' image_type]), image_type);

    close all
end

end
